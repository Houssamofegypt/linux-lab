#!/bin/sh
# compute_metrics.sh - binary classification metrics from a predictions CSV
# and a ground-truth CSV.
#
# Usage: compute_metrics.sh PREDICTIONS.csv TRUTH.csv
#
#   PREDICTIONS.csv  header: Peptide,Probability,Predicted Label
#   TRUTH.csv        header: Mutated_Epitope,Validation_Result
#
# Rows are paired by position (1st data row with 1st data row, ...), and the
# peptide in each pair must match. Columns are found by header name
# (case-insensitive), so column order does not matter. Blank lines and
# Windows (CRLF) line endings are accepted.
#
#   truth       = Validation_Result (0/1)
#   prediction  = Predicted Label   (0/1)  -> TP, TN, FP, FN, precision, recall, F1
#   score       = Probability              -> ROC AUC, PR AUC
#
# ROC AUC  area under the ROC curve; tied scores count as half, the same as
#          sklearn.metrics.roc_auc_score.
# PR AUC   reported two ways:
#          avg prec   average precision, the same as
#                     sklearn.metrics.average_precision_score
#          trapezoid  trapezoidal area under the precision-recall curve, the
#                     same as sklearn.metrics.auc(recall, precision)
#
# A metric that is undefined (e.g. precision with no positive predictions,
# or ROC AUC when the truth has only one class) is printed as NA.

set -eu

# Use "." as the decimal separator in awk and sort whatever the user's locale.
LC_ALL=C
export LC_ALL

prog=$(basename "$0")

usage() {
    echo "Usage: $prog PREDICTIONS.csv TRUTH.csv" >&2
    echo "  PREDICTIONS.csv: Peptide,Probability,Predicted Label" >&2
    echo "  TRUTH.csv:       Mutated_Epitope,Validation_Result" >&2
    exit 2
}

[ $# -eq 2 ] || usage
pred_file=$1
truth_file=$2

for f in "$pred_file" "$truth_file"; do
    if [ ! -f "$f" ] || [ ! -r "$f" ]; then
        echo "$prog: error: cannot read file '$f'" >&2
        exit 1
    fi
    if [ ! -s "$f" ]; then
        echo "$prog: error: file '$f' is empty" >&2
        exit 1
    fi
done

tmp=$(mktemp "${TMPDIR:-/tmp}/compute_metrics.XXXXXX")
trap 'rm -f "$tmp"' EXIT
trap 'exit 1' HUP INT TERM

# ---------------------------------------------------------------------------
# Step 1: pair the rows of the two files, validate them, and write one line
# per sample: "probability,truth,prediction".
# ---------------------------------------------------------------------------
awk -F',' -v prog="$prog" -v truth_name="$truth_file" -v pred_name="$pred_file" '
function clean(s) {
    gsub(/\r/, "", s)
    gsub(/^[ \t"]+|[ \t"]+$/, "", s)
    return s
}
function die(msg) {
    printf "%s: error: %s\n", prog, msg > "/dev/stderr"
    failed = 1
    exit 1
}
function is_number(s) {
    return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/
}
function label(s, what, file, line) {
    if (!is_number(s) || (s + 0 != 0 && s + 0 != 1))
        die(sprintf("%s line %d: %s must be 0 or 1, got \"%s\"", file, line, what, s))
    return s + 0
}
function column(file, name,    i) {
    for (i = 1; i <= NF; i++)
        if (tolower(clean($i)) == name)
            return i
    die(sprintf("%s: no \"%s\" column in header: %s", file, name, clean($0)))
}

# Header line of each file.
FNR == 1 {
    fileno++
    if (fileno == 1) {
        t_id = column(truth_name, "mutated_epitope")
        t_y  = column(truth_name, "validation_result")
    } else {
        p_id   = column(pred_name, "peptide")
        p_prob = column(pred_name, "probability")
        p_yhat = column(pred_name, "predicted label")
    }
    next
}

# Skip blank lines.
{ blank = $0; gsub(/[ \t\r]/, "", blank); if (blank == "") next }

# Truth file: remember peptide and label for each data row.
fileno == 1 {
    nt++
    t_pep[nt]  = clean($t_id)
    t_lab[nt]  = label(clean($t_y), "Validation_Result", truth_name, FNR)
    t_line[nt] = FNR
    next
}

# Predictions file: pair with the truth row at the same position.
fileno == 2 {
    np++
    if (np > nt)
        die(sprintf("%s has more data rows than %s (%d)", pred_name, truth_name, nt))
    pep  = clean($p_id)
    prob = clean($p_prob)
    yhat = label(clean($p_yhat), "Predicted Label", pred_name, FNR)
    if (pep != t_pep[np])
        die(sprintf("row %d: peptide mismatch: \"%s\" (%s line %d) vs \"%s\" (%s line %d)",
                    np, pep, pred_name, FNR, t_pep[np], truth_name, t_line[np]))
    if (!is_number(prob))
        die(sprintf("%s line %d: Probability is not a number: \"%s\"", pred_name, FNR, prob))
    print prob "," t_lab[np] "," yhat
}

END {
    if (failed) exit 1
    if (fileno < 2) die("could not read both files")
    if (np != nt)
        die(sprintf("row count mismatch: %s has %d data rows, %s has %d",
                    pred_name, np, truth_name, nt))
    if (np == 0) die("no data rows")
}
' "$truth_file" "$pred_file" > "$tmp"

# ---------------------------------------------------------------------------
# Step 2: sort by probability (highest first) and compute the metrics.
# Samples with the same probability form one threshold step on the curves.
# ---------------------------------------------------------------------------
sort -t',' -k1,1gr "$tmp" | awk -F',' '
# Close the current group of tied scores: record cumulative TP/FP counts.
function close_group() {
    ng++
    cum_tp[ng] = tp
    cum_fp[ng] = fp
}
function show(name, defined, value) {
    if (defined) printf "%-20s %.6f\n", name, value
    else         printf "%-20s NA\n", name
}

{
    score = $1 + 0; y = $2 + 0; yhat = $3 + 0

    # Confusion matrix from the predicted labels.
    if      (y == 1 && yhat == 1) TP++
    else if (y == 0 && yhat == 0) TN++
    else if (y == 0 && yhat == 1) FP++
    else                          FN++

    # Cumulative counts at each distinct score, for the ROC and PR curves.
    if (NR > 1 && score != prev) close_group()
    prev = score
    if (y == 1) tp++; else fp++
}

END {
    close_group()
    P = tp; N = fp      # number of truly positive / negative samples

    # ROC AUC: trapezoids between consecutive (FP, TP) points, starting at (0, 0).
    roc = 0; tp0 = 0; fp0 = 0
    for (k = 1; k <= ng; k++) {
        roc += (cum_fp[k] - fp0) * (cum_tp[k] + tp0) / 2
        tp0 = cum_tp[k]; fp0 = cum_fp[k]
    }

    # PR AUC: average precision = sum over thresholds of
    # (recall step) * precision, and the trapezoidal area starting at
    # (recall 0, precision 1).
    ap = 0; trapz = 0; r0 = 0; p0 = 1; tp0 = 0
    for (k = 1; k <= ng; k++) {
        if (P == 0) break
        r = cum_tp[k] / P
        p = cum_tp[k] / (cum_tp[k] + cum_fp[k])
        ap    += (cum_tp[k] - tp0) / P * p
        trapz += (r - r0) * (p + p0) / 2
        r0 = r; p0 = p; tp0 = cum_tp[k]
    }

    precision_ok = (TP + FP) > 0
    recall_ok    = (TP + FN) > 0
    f1_ok        = (2 * TP + FP + FN) > 0

    printf "%-20s %d (positives %d, negatives %d)\n", "Samples", NR, P, N
    printf "%-20s %d\n", "TP", TP
    printf "%-20s %d\n", "TN", TN
    printf "%-20s %d\n", "FP", FP
    printf "%-20s %d\n", "FN", FN
    show("Precision",          precision_ok, precision_ok ? TP / (TP + FP) : 0)
    show("Recall",             recall_ok,    recall_ok    ? TP / (TP + FN) : 0)
    show("F1 score",           f1_ok,        f1_ok        ? 2 * TP / (2 * TP + FP + FN) : 0)
    show("ROC AUC",            P > 0 && N > 0, (P > 0 && N > 0) ? roc / (P * N) : 0)
    show("PR AUC (avg prec)",  P > 0, ap)
    show("PR AUC (trapezoid)", P > 0, trapz)
}
'
