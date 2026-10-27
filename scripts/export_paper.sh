#!/usr/bin/env bash
set -euo pipefail
paper_project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "$paper_project_dir/envs/paper_python.sh" "$paper_project_dir/scripts/export_paper.py" "$@"
