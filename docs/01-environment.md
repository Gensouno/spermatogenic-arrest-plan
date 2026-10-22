# 01 环境可运行性报告

> 历史环境快照，非当前安装说明。部分路径、沙箱假设和依赖清单已经过时。当前入口见 [README](../README.md) 与 [复现说明](REPRODUCIBILITY.md)。

> 目标：给出**这台电脑上确实能跑**的环境方案。
> 本报告中的每一条"能装/能跑"都来自本机实测，不是通用教程搬运。
> 所有原始日志在 `logs/`（`logs/drive.log` 为安装汇总，`logs/i_<pkg>.log` 为单包日志）。

---

## 1. 机器事实（实测）

| 项目 | 实测值 | 取证方式 |
|---|---|---|
| 机型 / 架构 | Apple M4 / `arm64` | `uname -m` |
| CPU 核心 | 10 | `sysctl hw.ncpu` |
| 内存 | **16 GB**（17,179,869,184 B） | `sysctl hw.memsize` |
| 磁盘可用 | ~80 GB（跑完整个环境后仍剩 72 GB） | `df -h .` |
| 操作系统 | macOS 26.6.2 | `sw_vers` |
| R | 4.6.1 (2026-06-24) | `R.version.string` |
| Bioconductor | 3.23（BiocManager 1.30.27） | `BiocManager::version()` |
| Python | 3.14.7（conda base 的 `/opt/miniconda3/bin/python3`） | `python3 -V` |

> ⚠️ **纠正项目笔记**：`notes/思路整理.md:40` 写"这台 256G 的机器"，与实测的 16 GB 不符。
> 该笔记后续关于"下 SRA fastq 不现实"的判断在 16 GB 下**只会更成立**，结论方向不变，但数字必须按 16 GB 走。

---

## 2. R 库：项目内 `.Rlib`（实测 171 个包目录，934 MB）

系统库 `/Library/Frameworks/R.framework/Versions/4.6/Resources/library` **本机实测不可写**（见 §4.1），
因此所有新增包都装进项目内 `.Rlib`，运行时用：

```r
.libPaths(c(".Rlib", .libPaths()))
```

**实测加载结果：21/21 全部成功**（`logs/loadtest.log`，21 个包同时加载，耗时 38.8 s）。

### 关键包与用途分组（版本为实测）

| 分组 | 包（实测版本） |
|---|---|
| 差异表达 | limma 3.68.5、DESeq2 1.52.0、**edgeR 4.10.6** |
| 去卷积 | **nnls 1.6**、**quadprog 1.5.8**、EpiDISH 2.28.0、TOAST 1.26.0 |
| 富集 / 打分 | clusterProfiler 4.20.0、**fgsea 1.38.0**、**GSVA 2.6.12**、**topGO 2.64.0** |
| 可视化 | pheatmap 1.0.13、**ComplexHeatmap 2.28.0**、UpSetR 1.4.1、corrplot 0.95、ggplot2 4.0.3、patchwork 1.3.2 |
| 单细胞 | Seurat 5.6.0、SingleCellExperiment 1.34.0、SingleR 2.14.2、scran 1.40.0、scater 1.40.2 |
| 其它 | **sva 3.60.0**、**WGCNA 1.74**、variancePartition 1.42.0、muscat 1.26.0、speckle 1.12.0、lme4 2.0.6、**openxlsx 4.2.9** |

（加粗 = 本次补齐装入 `.Rlib` 的包。）

### 2.1 逐个包实测结果表（本项目需要的包）

`耗时` = 该包及其依赖的**完整安装**耗时（已存在的包记为 SKIP）。`RSS` = 安装进程最大常驻内存。
`来源` = 实测命中的下载仓库；**全部为 arm64 预编译二进制，无一编译**（见 §4.7）。

| 包名 | 来源 | 主要依赖（实测） | 耗时 | RSS | 状态 |
|---|---|---|---|---|---|
| edgeR | Bioconductor 3.23 | 系统已有 limma/locfit | 35.9 s | 247 MiB | ✅ 装入 `.Rlib` |
| sva | Bioconductor 3.23 | xtable, annotate, genefilter | 39.9 s | 249 MiB | ✅ |
| fgsea | Bioconductor 3.23 | cowplot, fastmatch | 44.9 s | 238 MiB | ✅ |
| GSVA | Bioconductor 3.23 | 19 个：HDF5Array, GSEABase, SpatialExperiment, magick… | 59.4 s | 235 MiB | ✅ |
| ComplexHeatmap | Bioconductor 3.23 | rjson, circlize, GetoptLong, colorspace, clue, GlobalOptions, doParallel | 41.1 s | 240 MiB | ✅ |
| nnls | CRAN（arm64 二进制） | 无 | 未测† | 未测† | ✅ **会话开始前已存在于 `.Rlib`** |
| quadprog | CRAN（arm64 二进制） | 无 | 2.6 s* | 230.5 MiB*（三包合计） | ✅ 首批装入 |
| preprocessCore | Bioconductor 3.23 | 无 | 35.1 s | 234.3 MiB | ✅ |
| openxlsx | CRAN（arm64 二进制） | 无 | 8.8 s* | 同上（合计） | ✅ 首批装入 |
| corrplot | CRAN（arm64 二进制） | 无 | 14.0 s* | 同上（合计） | ✅ 首批装入 |
| SingleCellExperiment | Bioconductor 3.23 | SummarizedExperiment 链 | — | 85 MiB | ✅（scran/GSVA 依赖带装） |
| SingleR | Bioconductor 3.23 | beachmat, assorthead | 36.1 s | 247 MiB | ✅ |
| TOAST | Bioconductor 3.23 | proxy, forcats, e1071, locfdr, ggstats, EpiDISH, corpcor, GGally | 37.0 s | 248 MiB | ✅ |
| WGCNA | CRAN | backports, checkmate, htmlTable, Formula, dynamicTreeCut, fastcluster, Hmisc, impute | 50.1 s | 247 MiB | ✅ |
| variancePartition | Bioconductor 3.23 | 31 个：lme4, lmerTest, pbkrtest, RcppArmadillo, forecast, remaCor… | 46.4 s | 237 MiB | ✅ |
| muscat | Bioconductor 3.23 | 27 个：glmmTMB, TMB, RTMB, scater, blme, BiocSingular… | 64.7 s | 236 MiB | ✅ |
| speckle | Bioconductor 3.23 | 49 个：Seurat, spatstat.*, plotly, shiny, future…（**依赖最重**） | 56.0 s | 250 MiB | ✅ |
| scran | Bioconductor 3.23 | bluster, metapod | 44.1 s | 238 MiB | ✅ |
| scater | Bioconductor 3.23 | 随 muscat 带装 | — | 85 MiB | ✅ |
| topGO | Bioconductor 3.23 | 系统已有 graph | 49.3 s | 250 MiB | ✅ |
| remotes | CRAN（arm64 二进制） | 无 | 11.1 s | 198 MiB | ✅（供后续 GitHub 装包，但本次 GitHub 不可达） |

\* 该包在较早的**首批验证批次**中已装入 `.Rlib`，故后续驱动脚本将其记为 `SKIP(exists)`；耗时/RSS 取自首次真实安装。
† `nnls` 在会话开始前就已存在于 `.Rlib`，本次未重新安装，故无实测耗时/RSS。

**失败项**：

| 包名 | 来源 | 耗时 | 状态 | 替代方案 |
|---|---|---|---|---|
| MuSiC | GitHub `xuranw/MuSiC` | 6 次重试全败（约 40 min） | ❌ 下载截断 / DNS 失败 | `nnls` + `quadprog` 自行实现 NNLS |
| bisque | GitHub `cozygene/bisque` | 6 次重试全败 | ❌ `Could not resolve host: codeload.github.com` | `EpiDISH`（已装） |
| BisqueRNA | 不在 CRAN/Bioc | — | ❌ 分发渠道不存在 | 同上 |
| granulator | 不在 CRAN/Bioc | — | ❌ 分发渠道不存在 | 不纳入方案 |
| glmTMB | — | — | — | **该包名不存在**；实际需要的是 `glmmTMB`（已随 muscat 装入） |

---

## 3. Python 环境：项目内 `.venv`

```
.venv/bin/python  →  Python 3.14.7
numpy 2.5.4 | pandas 3.0.6 | scipy 1.18.1
scikit-learn 1.9.1 | openpyxl 3.1.5 | matplotlib 3.11.2
```

**实测能读写的矩阵格式**（写入 → 读回 → 校验形状，`logs/setup_full.log`）：

| 格式 | 结果 |
|---|---|
| `.xlsx`（openpyxl 引擎） | ✅ 正常 |
| `.csv.gz`（gzip 压缩） | ✅ 正常 |
| `.txt.gz`（tab 分隔 + gzip） | ✅ 正常 |
| `.tsv` 纯文本 | ✅ 正常（用真实文件 `data/geo_candidates.tsv`，448×14） |

---

## 4. 本机做不到的事 + 替代路径

### 4.1 系统 R library 不可写 —— ✅ 结论成立，但**原因需要修正**

**实测错误**：
```
Warning: cannot create file
'/Library/Frameworks/R.framework/Versions/4.6/Resources/library/__wtest__',
reason 'Operation not permitted'
```

**替代**：项目内 `.Rlib` + `.libPaths(c(".Rlib", .libPaths()))`。

> **精确化（与"已实测约束"的差异）**：该目录实测属性为
> `owner=root group=admin mode=drwxrwxr-x`，而当前用户的组为 `admin`（`id` 可见 `80(admin)`）。
> 也就是说**按 Unix 权限，组写位本应允许写入**；实测返回的是 `EPERM`（errno=1）而非 `EACCES`（13）。
> 判定：**在本机当前配置下该目录不可写**（这是运行 R 脚本的实际场景，方案按此设计），
> 但根因是 **本机的写权限限制**（当前配置只允许写 `$HOME/Downloads`），而非 macOS 权限本身。
> 若将来在 此受限环境之外的普通终端里运行，该目录**可能**是可写的 —— 但**不要依赖这一点**，
> `.Rlib` 方案在任何情况下都正确，且能保证环境可复现。

### 4.2 conda —— ❌ **结论需要推翻：conda 并非"完全不可用"**

两个报错**都能复现**（原文照抄）：

```
NoWritableEnvsDirError: No writeable envs directories configured.
  - $HOME/.conda/envs
  - /opt/miniconda3/envs
```
```
CondaToSPermissionError: Unable to read/write path
($HOME/Library/Caches/conda-anaconda-tos/3c9d068aa053e2a1c4313fe3391b7a8ee57c4fbd09c4e8aeae49ef333a740150.cache)
for https://repo.anaconda.com/pkgs/main. Please check permissions.
```

**但是**：这两个错误都是**默认路径**造成的，把路径重定向进项目后 conda 可以正常工作，已实测建成可用环境：

```bash
export XDG_CACHE_HOME="$PWD/.cache/xdg"        # 解决 CondaToSPermissionError
export CONDA_PKGS_DIRS="$PWD/.cache/conda-pkgs" # 解决 NoWritablePkgsDirError
export CONDA_REGISTER_ENVS=false                # 解决 ~/.conda/environments.txt 不可写
conda create -p "$PWD/.cache/conda-envs/probe" -c conda-forge --override-channels python=3.12 -y
```

**实测结果（A/B 对照）**：

| 条件 | 命令 | 结果 |
|---|---|---|
| **A：带重定向** | 上述 3 个 export + `conda create -p .cache/conda-verify -c conda-forge --override-channels python=3.12 -y` | `EXIT=0`，`.cache/conda-verify/bin/python -V` → **Python 3.12.15** ✅，`import sys` 正常 |
| **B：不带重定向（对照）** | `env -u XDG_CACHE_HOME -u CONDA_PKGS_DIRS -u CONDA_REGISTER_ENVS conda create ...` | `CondaToSPermissionError` ❌ |

对照实验证实：**失败完全由默认缓存路径引起，与 conda 本身无关**。
（注：A 组验证用的临时前缀已在验证后清理，命令可原样重跑复现。）

> **证据链**：不带 `CONDA_REGISTER_ENVS=false` 时，conda 把事务**完整执行完**，
> 最后只在写 `~/.conda/environments.txt` 时失败并回滚（日志：`Executing transaction: ...done` → `Rolling back transaction: ...done`
> → `[Errno 1] Operation not permitted: '$HOME/.conda/environments.txt'`）。
> 这说明失败点只是一个注册文件，**不是 conda 本身不可用**。

**结论与建议**：
- 早期判断"conda 完全不可用"**不准确**。
- **但方案的最终选择仍然是"不用 conda"** —— 理由从"conda 不可用"改为：
  1. 项目已有可用的 `.venv` + `.Rlib`，切换成本为零；
  2. 受限环境下 conda 需要对 4 个环境变量做非默认重定向，脆弱且难排查；
  3. 本项目不需要 conda 的二进制分发能力（R 包已全部拿到 arm64 预编译二进制）。
- 该 workaround 已记录在此，**供在此受限环境之外的普通终端里自行取用**；`envs/setup.sh` 中仍以 venv 为主路径。

### 4.3 matplotlib 配置目录不可写 —— ✅ 结论成立

**实测警告**：
```
$HOME/.matplotlib is not a writable directory
Matplotlib created a temporary cache directory at /var/folders/.../matplotlib-g2pkr2gc ...
```

**替代**：`export MPLCONFIGDIR="$PWD/.mplcache"`（`.mplcache/` 已建；已写入 `envs/setup.sh` 的 `MPLCONFIGDIR`）。

### 4.4 pip 缓存不可写 —— ✅ 结论成立（可忽略）

**实测警告**：
```
WARNING: The directory '$HOME/Library/Caches/pip' or its parent directory is not owned
or is not writable by the current user. The cache has been disabled.
```
**替代**：`export PIP_CACHE_DIR="$PWD/.cache/pip"`（`setup.sh` 已设）。不设也只是每次多一条 warning，不影响安装。

### 4.5 GitHub 装包不可靠 —— ✅ 结论成立（实测失败）

`MuSiC`（`xuranw/MuSiC`）与 `bisque`（`cozygene/bisque`）实测**全部失败**，6 次重试无一成功：

| 尝试 | 结果 |
|---|---|
| 1 | `curl: (92) HTTP/2 stream 1 was not closed cleanly: CANCEL` → 18,621,917 B 截断 |
| 2 | `curl: (92) ... CANCEL` → 14,812,891 B 截断 |
| 3 | `curl: (28) Operation too slow. Less than 500 bytes/sec ... last 120 seconds` → 16,905,763 B 截断 |
| 4–6 | `curl: (6) Could not resolve host: codeload.github.com` |
| bisque 1–6 | `curl: (6) Could not resolve host: codeload.github.com` |

- `gzip -t` 判定下载物 **INVALID/TRUNCATED**；`git clone --depth 1` 同样失败
  （`error: RPC failed; curl 18 transfer closed with outstanding read data remaining` / `fatal: early EOF`）。
- 残余文件 `envs/vendor/MuSiC.tar.gz`（16.9 MB）**是坏文件，不要使用**。

**替代**：去卷积用 CRAN 的 **`nnls` 1.6 + `quadprog` 1.5.8** 自行实现 NNLS/约束最小二乘，
配合已装好的 **EpiDISH 2.28.0**（CRAN/Bioc 分发，实测可用）作为参考实现。**不依赖 GitHub-only 的 MuSiC / BisqueRNA**。

### 4.6 网络速率 —— ⚠️ **与早期估算不符，以实测为准**

早期估算"CRAN/GitHub 可达但慢""CRAN 索引可拉但下载慢（实测 192 B/s）"。**实测远好于此**：

| 目标 | 实测 |
|---|---|
| CRAN 索引 `PACKAGES.gz` | 1,945,195 B / 6.30 s = **309 KB/s** |
| CRAN 包（源码 `corrplot_0.95.tar.gz`） | **1,015,708 B/s** |
| CRAN 包（源码 `quadprog_1.5-8.tar.gz`） | 39 KB/s |
| **CRAN 实际使用的 arm64 二进制**（`bin/macosx/sonoma-arm64/contrib/4.6/corrplot_0.95.tgz`） | **1,545,603 B/s** |
| PyPI 直连 | `curl: (56) Recv failure: Connection reset by peer`（**不可用**） |
| 清华 PyPI 镜像 | **2,647,987 B/s** |
| Bioconductor 索引 | 236,625 B / 3.53 s = 67 KB/s |

**结论**：CRAN 与 PyPI-清华镜像**速度良好**，不是瓶颈。真正的瓶颈是 **GitHub**（§4.5）与 **NCBI GEO**
（早期估算的"单流 16–35 KB/s、8 并发最优"本次仅部分复现：单流实测 14 KB/s，与区间吻合；
并发标定因超时中断，未取得可靠数据，**方案按 8 并发**沿用早期结论，但标注为"未独立复核"）。

### 4.7 编译 —— ❌ **结论需要修正：实际上根本没编译**

早期版本称"R 包源码编译可行"、`setup.sh` 初版也写了"Bioconductor 走源码编译"。
**实测：21 个包没有任何一个发生编译**。

- 全部日志中**没有任何**编译痕迹（`grep -lE "clang|\.cpp|\*\* libs"` → 空）。
- CRAN 走的是 `bin/macosx/sonoma-arm64/contrib/4.6/` **预编译二进制**。
- Bioconductor 3.23 **同样提供 arm64 二进制**：
  `https://bioconductor.org/packages/3.23/bioc/bin/macosx/sonoma-arm64/contrib/4.6/PACKAGES.gz` → HTTP 200（227,282 B）。
- `edgeR` 的 `DESCRIPTION` 实测：`Built: R 4.6.1; aarch64-apple-darwin23`；目录内只有 `libs/edgeR.so`，无 `.c/.cpp`。

**影响**：安装速度远快于预期（单包 35–65 s，非"编译几十分钟"），且**不依赖 clang/CommandLineTools**。
`setup.sh` 中的相关注释已同步更正。

> 注：`glmTMB` 这个包**不存在**；`muscat` 实际拉取的是 **`glmmTMB`**（实测已随 muscat 装入）。

---

## 5. 实测性能与 16 GB 内存预算

### 5.1 实测内存峰值（`/usr/bin/time -l` 的 `maximum resident set size`，MiB = 2²⁰ B）

| 场景 | 原始字节 | 实测峰值 |
|---|---|---|
| 小包（nnls/quadprog/openxlsx/corrplot）安装 | 88,064,000 | **84 MiB** |
| 单个 R 包安装（含依赖，R 进程基线） | ~247,000,000 | **234–250 MiB** |
| edgeR 20000 基因 × 40 样本 QL 检验 | 228,982,784 | **218 MiB** |
| WGCNA + SingleR + scran/TOAST 等组合加载 | 1,091,125,248 | **1,041 MiB** |
| fgsea + ComplexHeatmap + sva + GSVA 组合 | 1,220,214,784 | **1,164 MiB** |
| `setup.sh` 完整跑一遍 | 1,460,502,528 / 1,492,975,616 | **1,393 / 1,424 MiB** |
| 21 个目标包**同时** `library()` 加载 | 1,501,085,696 | **1,432 MiB** |

### 5.2 分阶段内存上限估算（16 GB = 16,384 MiB 预算）

| 阶段 | 实测/估算峰值 | 占 16 GB | 判断 |
|---|---|---|---|
| 环境安装（setup.sh） | 1,393 MiB | 8.5% | 安全，可与其它任务并行 |
| 单数据集 bulk DE（20k × 40） | 218 MiB | 1.3% | 安全 |
| 多数据集合并 + sva 校正（~20k × 200） | ~1,000–2,000 MiB | 6–12% | 安全 |
| 富集分析（fgsea/GSVA，200 基因集） | 1,164 MiB | 7.1% | 安全 |
| WGCNA（3k–20k 基因） | 1,041 MiB – ~4 GB | 6–25% | 基因数 >15k 时需注意 |
| **单细胞去卷积 / SingleR（>5 万细胞）** | **4–8 GB** | **25–50%** | ⚠️ 主要风险点 |
| **Seurat 全量单细胞对象（>10 万细胞）** | **>8 GB** | **>50%** | ⚠️ 需分块或降采样 |

**结论**：bulk RNA-seq 全流程在 16 GB 下**舒适**（峰值 <15%）。
**唯一需要管控的是单细胞阶段** —— 建议：
1. 用 `SingleCellExperiment` + `HDF5Array`（已随 GSVA 装入）做**磁盘后端**，避免全量入内存；
2. `SingleR` 参考集分块比对；
3. 超过 10 万细胞时**先降采样**再估计细胞比例，或在 8 并发下载后**逐数据集流式处理**，不合并成单个大对象。

---

## 6. 一键安装

```bash
cd "spermatogenic-arrest-plan"
bash envs/setup.sh          # 完整安装（幂等，可重复执行）
bash envs/setup.sh --check  # 只体检，不安装
```

**实测幂等性**：重复运行两次（所有包已存在）分别耗时 **14.9 s** 与 **10.1 s**、均 `SETUP-EXIT=0`，
全部包走 `SKIP(exists)` 分支，验证段仍输出 `R 包加载: 21/21 成功`。

脚本会自动完成：项目内缓存重定向（`XDG_CACHE_HOME` / `CONDA_PKGS_DIRS` / `CONDA_REGISTER_ENVS=false` /
`MPLCONFIGDIR` / `PIP_CACHE_DIR`）→ 硬件自检 → `.Rlib` 装包 → `.venv` 装包 → **加载验证** → **读写验证**。

---

## 7. 可删除的残留（仅列出，**未删除**）

| 路径 | 大小 | 说明 |
|---|---|---|
| `.venvtest/` | 47 MB | 废弃的测试环境，只有 numpy |
| `.conda_meta/` | 151 MB | conda workaround 试验残留（验证用的 probe 前缀已清理） |
| `.condapkgs/` | 15 MB | 早先的 conda 包缓存 |
| `.condatmp/` | 8 KB | conda 临时目录 |
| `envs/vendor/` | 17 MB | **损坏的** MuSiC tarball（§4.5），无保留价值 |

合计约 **230 MB**，删除不影响 `setup.sh` 的可复现性。

> **注意：`.cache/`（151 MB）不是残留，请保留。** 它是 `setup.sh` 主动使用的项目内缓存
> （`XDG_CACHE_HOME` / `CONDA_PKGS_DIRS` / `PIP_CACHE_DIR`），删掉只会让下次运行重新下载。

---

## 8. 与"已实测约束"的差异汇总

| # | 初始设计约束 | 实测结论 | 是否影响方案 |
|---|---|---|---|
| 1 | conda 完全不可用 | **不准确**：重定向 4 个环境变量后 conda 可用（已建成 Python 3.12.15 环境）。方案**仍不用 conda**，但理由改为工程取舍 | 否（结论不变，理由更正） |
| 2 | 系统 R library 不可写 | **在本机当前配置下成立**；根因是 本机写权限限制而非 Unix 权限（当前用户 ∈ `admin`，目录组写位为 `rwx`） | 否（`.Rlib` 方案照旧） |
| 3 | 当前环境只允许写 `$HOME/Downloads` | **成立**。补充：`/tmp` 也可写；`~` 下不可写 | 否 |
| 4 | R 包需源码编译 | **不成立**：CRAN 与 Bioconductor 均命中 arm64 预编译二进制，**零编译** | 是（安装更快，不依赖 clang） |
| 5 | CRAN 下载慢（192 B/s） | **不成立**：索引 309 KB/s、二进制 1.5 MB/s | 是（安装时间预期大幅下调） |
| 6 | PyPI 可用 | **部分不成立**：PyPI 直连被重置，须用清华镜像（2.6 MB/s） | 是（`setup.sh` 已内置镜像） |
| 7 | GitHub 装 MuSiC/Bisque | **成立**：实测 6 次全败，DNS/HTTP2 双重故障 | 否（改用 nnls/quadprog） |
| 8 | 内存 16 GB（笔记写 256G） | **16 GB 成立**，笔记有误 | 是（单细胞阶段需管控） |
| 9 | GEO 8 并发最优 | **未独立复核**（标定超时中断），沿用原结论 | 否 |
