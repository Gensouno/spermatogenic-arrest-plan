#!/usr/bin/env bash
# Minimal audited pipeline dependencies; never moves Python/R or uses sudo.
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"
command -v Rscript >/dev/null || { echo 'Rscript not found; install R before setting up this project.' >&2; exit 1; }
if [ "${1:-}" != "--check" ]; then
  if [ ! -x .venv/bin/python ]; then
    "${PYTHON_BIN:-python3}" -c 'import sys; sys.exit(0 if sys.version_info >= (3,12) else "Python >=3.12 required; select an interpreter with PYTHON_BIN")'
    "${PYTHON_BIN:-python3}" -m venv .venv
  fi
  .venv/bin/python -c 'import sys; sys.exit(0 if sys.version_info >= (3,12) else "Existing environment needs Python >=3.12; rebuild explicitly, not in place")'
  bash envs/python.sh -m pip install -r envs/requirements-analysis.txt
  mkdir -p .Rlib
  bash envs/rscript.sh -e '
    lib <- file.path(getwd(), ".Rlib")
    options(repos=c(CRAN="https://cloud.r-project.org"))
    if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",lib=lib)
    pkgs <- c("data.table","edgeR","limma","AnnotationDbi","org.Mm.eg.db",
              "affy","mouse4302.db","mouse4302cdf","statmod")
    missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
    if(length(missing)) BiocManager::install(missing,lib=lib,ask=FALSE,update=FALSE)
  '
fi
bash envs/python.sh -c '
from importlib.metadata import version
import sys
if sys.version_info < (3, 12): raise SystemExit("Python >=3.12 required for the tested package set")
print(sys.version)
for package in ["numpy", "pandas", "scipy", "matplotlib", "openpyxl"]:
    __import__(package)
    print(package, version(package))
'
bash envs/rscript.sh -e '
  if(getRversion()<"4.4") stop("R >=4.4 required; package versions depend on matching Bioconductor release")
  pkgs <- c("data.table","edgeR","limma","AnnotationDbi","org.Mm.eg.db",
            "affy","mouse4302.db","mouse4302cdf","statmod")
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Missing R packages: ",paste(missing,collapse=", "))
  print(sapply(pkgs,function(p)as.character(packageVersion(p))))
  print(R.version.string)
'
