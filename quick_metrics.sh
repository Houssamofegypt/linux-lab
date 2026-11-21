#!/bin/sh
# quick_metrics.sh PRED.csv TRUTH.csv  -- quick & dirty, minimal checks.
# PRED cols (fixed order):  Peptide,Probability,Predicted Label
# TRUTH cols (fixed order): Mutated_Epitope,Validation_Result
# Rows are paired by line position. PR_AUC = average precision.
export LC_ALL=C
paste -d, "$1" "$2" | tr -d '\r' | tail -n +2 | sort -t, -k2,2gr | awk -F, '
{ s = $2; yh = $3; y = $5
  if ($1 != $4) bad++
  if (y && yh) TP++; else if (!y && !yh) TN++; else if (yh) FP++; else FN++
  if (NR > 1 && s != prev) { n++; T[n] = tp; F[n] = fp }
  prev = s; if (y) tp++; else fp++ }
END { n++; T[n] = tp; F[n] = fp; P = tp; N = fp
  for (k = 1; k <= n; k++) {
    roc += (F[k] - F[k-1]) * (T[k] + T[k-1]) / 2
    ap  += (T[k] - T[k-1]) / P * T[k] / (T[k] + F[k]) }
  if (bad) printf "WARNING: %d rows with mismatched peptides\n", bad
  printf "TP=%d TN=%d FP=%d FN=%d\n", TP, TN, FP, FN
  printf "precision=%.4f recall=%.4f F1=%.4f\n", TP/(TP+FP), TP/(TP+FN), 2*TP/(2*TP+FP+FN)
  printf "ROC_AUC=%.4f PR_AUC=%.4f\n", roc/(P*N), ap }'
