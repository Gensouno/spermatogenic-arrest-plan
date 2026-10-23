#!/usr/bin/env bash
# Sourced by project launchers; no machine-specific proxy is stored in Git.
if [ -f "$PROJECT_DIR/.download-proxy.local" ]; then
  project_download_proxy="$(< "$PROJECT_DIR/.download-proxy.local")"
  case "$project_download_proxy" in
    http://?*|https://?*) ;;
    *) echo 'Invalid .download-proxy.local: expected one HTTP(S) proxy URL.' >&2; exit 1 ;;
  esac
  if [[ "$project_download_proxy" =~ [[:space:]] ]]; then
    echo 'Invalid .download-proxy.local: whitespace is not allowed.' >&2
    exit 1
  fi
  export HTTP_PROXY="$project_download_proxy" HTTPS_PROXY="$project_download_proxy"
  export http_proxy="$project_download_proxy" https_proxy="$project_download_proxy"
  export ALL_PROXY="$project_download_proxy" all_proxy="$project_download_proxy"
  export NO_PROXY='localhost,127.0.0.1,::1' no_proxy='localhost,127.0.0.1,::1'
  unset project_download_proxy
fi
