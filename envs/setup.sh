#!/usr/bin/env bash
# =============================================================================
#  生精阻滞共性分析 —— 项目本地环境搭建（幂等，可重复执行）
# =============================================================================
#  实测机器: Apple M4, 10 核, 16 GB RAM, macOS (darwin arm64)
#  可用磁盘: ~81 GB      R 4.6.1 + Bioconductor 3.23
#
#  本脚本的三个硬性前提（都已实测确认，见 docs/01-environment.md）:
#    1. 系统 R library 不可写 —— 所有 R 包装进项目内 .Rlib
#    2. conda 默认路径不可用 —— 必须把 envs/pkgs/cache 全部重定向进项目
#    3. 当前环境只允许写项目目录 —— 所有缓存与临时目录都指向项目内
#
#  用法:  bash envs/setup.sh          # 在项目根目录执行
#         bash envs/setup.sh --check  # 只检查、不安装
# =============================================================================
set -uo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ" || exit 1
CHECK_ONLY=0
[ "${1:-}" = "--check" ] && CHECK_ONLY=1

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }
ok()  { printf '   [OK]   %s\n' "$*"; }
bad() { printf '   [FAIL] %s\n' "$*"; }

mkdir -p logs envs reports .Rlib
mkdir -p data/tmp

# =============================================================================
#  0. 环境变量：把每一个缓存目录都钉死在项目内
# =============================================================================
# 为什么必须这么做（实测）:
#   - conda 默认找 $HOME/.conda/environments.txt 与
#     $HOME/Library/Caches/conda-anaconda-tos/ ；这两处都在可写目录之外，
#     任何写入都返回 EPERM(errno=1)，conda 会在"事务已执行完"之后回滚。
#   - pip/matplotlib 同理：~/Library/Caches、~/.matplotlib 都不可写。
# 把这些路径重定向到项目内，conda 与 pip 都能正常工作。
export XDG_CACHE_HOME="$PROJ/.cache/xdg"
export CONDA_PKGS_DIRS="$PROJ/.cache/conda-pkgs"
export CONDA_ENVS_DIRS="$PROJ/.cache/conda-envs"
export CONDA_REGISTER_ENVS=false      # 关键：不写 ~/.conda/environments.txt
export MPLCONFIGDIR="$PROJ/.mplcache"   # 实测 ~/.matplotlib 不可写，必须重定向
export PIP_CACHE_DIR="$PROJ/.cache/pip"
mkdir -p "$XDG_CACHE_HOME" "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS" "$MPLCONFIGDIR" "$PIP_CACHE_DIR"

# 国内镜像（实测 PyPI 直连被重置，清华镜像 2.6 MB/s）
PIP_INDEX="https://pypi.tuna.tsinghua.edu.cn/simple"

# =============================================================================
#  1. 硬件与前置条件自检
# =============================================================================
say "1/6 硬件与前置条件"
MEM_BYTES=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
MEM_GB=$(( MEM_BYTES / 1073741824 ))
echo "   内存: ${MEM_GB} GB   核心: $(sysctl -n hw.ncpu 2>/dev/null)   磁盘可用: $(df -h "$PROJ" | tail -1 | awk '{print $4}')"
if [ "$MEM_GB" -lt 16 ]; then
  bad "内存 ${MEM_GB} GB 低于方案的 16 GB 假设"
else
  ok "内存 ${MEM_GB} GB"
fi
for t in R Rscript python3; do
  command -v "$t" >/dev/null 2>&1 && ok "$t -> $(command -v "$t")" || bad "$t 未找到"
done
Rscript -e 'cat("   R:", R.version.string, "| Bioc:", as.character(BiocManager::version()), "\n")' 2>/dev/null

# =============================================================================
#  2. R 本地库 .Rlib
# =============================================================================
say "2/6 R 本地库 .Rlib"
mkdir -p "$PROJ/.Rlib"
# 验证系统库确实不可写（不是猜测）
Rscript -e '
  sys <- .libPaths()[length(.libPaths())]
  r <- tryCatch({ f <- file.path(sys, ".wtest"); file.create(f); file.remove(f); "WRITABLE" },
                error=function(e) "NOT-WRITABLE", warning=function(w) "NOT-WRITABLE")
  cat("   system lib:", sys, "->", r, "\n")
' 2>/dev/null
echo "   项目本地库: $PROJ/.Rlib"

# 需要装入 .Rlib 的包（实测 21/21 全部加载成功）
R_PKGS="edgeR sva fgsea GSVA ComplexHeatmap UpSetR nnls quadprog preprocessCore \
WGCNA SingleCellExperiment SingleR TOAST variancePartition muscat speckle scran \
scater topGO corrplot openxlsx remotes"

if [ "$CHECK_ONLY" -eq 0 ]; then
  say "3/6 安装 R 包（串行；CRAN 与 Bioconductor 均命中 arm64 预编译二进制，实测无需编译）"
  for p in $R_PKGS; do
    Rscript envs/install_pkgs.R "$p" 2>&1 | grep -E '^\[RESULT\]|^\[ERROR\]' || echo "   [WARN] $p 无结果行"
  done
else
  say "3/6 安装 R 包（--check，跳过）"
fi

# =============================================================================
#  3. Python .venv
# =============================================================================
say "4/6 Python 虚拟环境 .venv"
if [ ! -x "$PROJ/.venv/bin/python" ]; then
  [ "$CHECK_ONLY" -eq 0 ] && python3 -m venv "$PROJ/.venv" && ok "已创建 .venv"
else
  ok ".venv 已存在: $("$PROJ/.venv/bin/python" -V 2>&1)"
fi
if [ "$CHECK_ONLY" -eq 0 ]; then
  # 实测: PyPI 直连 connection reset；清华镜像正常
  "$PROJ/.venv/bin/python" -m pip install --quiet -i "$PIP_INDEX" \
      numpy pandas scipy scikit-learn openpyxl matplotlib \
    && ok "Python 依赖已安装" || bad "Python 依赖安装失败"
fi

# =============================================================================
#  4. 验证（这一步是"确实能跑"的证据，不是装饰）
# =============================================================================
say "5/6 验证 R 包可加载"
Rscript -e '
.libPaths(c(".Rlib", .libPaths()))
pkgs <- strsplit(trimws("edgeR sva fgsea GSVA ComplexHeatmap UpSetR nnls quadprog
preprocessCore WGCNA SingleCellExperiment SingleR TOAST variancePartition muscat
speckle scran scater topGO corrplot openxlsx"), "[[:space:]]+")[[1]]
ok <- 0; bad <- character()
for (p in pkgs) {
  r <- tryCatch({ suppressMessages(library(p, character.only=TRUE)); TRUE },
                error=function(e) FALSE)
  if (r) ok <- ok + 1 else bad <- c(bad, p)
}
cat(sprintf("   R 包加载: %d/%d 成功\n", ok, length(pkgs)))
if (length(bad)) cat("   失败:", paste(bad, collapse=", "), "\n")
quit(status = if (ok == length(pkgs)) 0 else 1)
' 2>&1

say "6/6 验证 Python 能读 xlsx / csv.gz / txt.gz"
"$PROJ/.venv/bin/python" - <<'PY' 2>&1
import os, sys
os.environ.setdefault("MPLCONFIGDIR", os.path.abspath(".cache/matplotlib"))
import pandas as pd, numpy as np
df = pd.DataFrame({"gene": ["G1","G2"], "v": [1.0, 2.0]})
os.makedirs("data/tmp", exist_ok=True)
df.to_excel("data/tmp/_probe.xlsx", index=False); assert pd.read_excel("data/tmp/_probe.xlsx").shape == (2,2)
df.to_csv("data/tmp/_probe.csv.gz", index=False, compression="gzip"); assert pd.read_csv("data/tmp/_probe.csv.gz").shape == (2,2)
df.to_csv("data/tmp/_probe.txt.gz", sep="\t", index=False, compression="gzip"); assert pd.read_csv("data/tmp/_probe.txt.gz", sep="\t").shape == (2,2)
for f in ("_probe.xlsx","_probe.csv.gz","_probe.txt.gz"): os.remove(f"data/tmp/{f}")
print("   xlsx / csv.gz / txt.gz 读写均正常")
import numpy, scipy, sklearn, matplotlib, openpyxl
print("   numpy", numpy.__version__, "| pandas", pd.__version__, "| scipy", scipy.__version__)
print("   sklearn", sklearn.__version__, "| openpyxl", openpyxl.__version__, "| matplotlib", matplotlib.__version__)
PY

say "完成"
echo "   R 本地库 : $PROJ/.Rlib"
echo "   Python   : $PROJ/.venv"
echo "   缓存      : $PROJ/.cache"
echo
echo "   日常使用（每个 R 脚本开头都要加这一行）:"
echo "     .libPaths(c(\"$PROJ/.Rlib\", .libPaths()))"
