#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$PROJECT_DIR/.paper-tools/packages${PYTHONPATH:+:$PYTHONPATH}"
exec bash "$PROJECT_DIR/envs/python.sh" "$@"
