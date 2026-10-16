#!/usr/bin/env bash
cd "$(dirname "$0")/.."
PKGS="ComplexHeatmap preprocessCore openxlsx corrplot nnls quadprog GSVA WGCNA SingleCellExperiment SingleR TOAST variancePartition muscat speckle scran scater topGO"
for p in $PKGS; do
  echo "##### $p  $(date +%H:%M:%S)"
  /usr/bin/time -l Rscript envs/install_pkgs.R "$p" > "logs/i_$p.log" 2>&1
  grep -E "^\[(RESULT|ERROR)\]" "logs/i_$p.log" || echo "[RESULT] $p NO-LINE"
  grep -E "maximum resident set size" "logs/i_$p.log" | sed 's/^/    rss:/'
done
echo "##### DRIVE-DONE $(date +%H:%M:%S)"
