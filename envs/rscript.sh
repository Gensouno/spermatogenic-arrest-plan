#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"
source "$PROJECT_DIR/envs/network.sh"
export R_LIBS_USER="$PROJECT_DIR/.Rlib${R_LIBS_USER:+:$R_LIBS_USER}"
exec Rscript "$@"
