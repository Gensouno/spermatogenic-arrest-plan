#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"
source "$PROJECT_DIR/envs/network.sh"
export MPLCONFIGDIR="$PROJECT_DIR/.mplcache"
export XDG_CACHE_HOME="$PROJECT_DIR/.cache/xdg"
export PIP_CACHE_DIR="$PROJECT_DIR/.cache/pip"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$PIP_CACHE_DIR"
exec "$PROJECT_DIR/.venv/bin/python" "$@"
