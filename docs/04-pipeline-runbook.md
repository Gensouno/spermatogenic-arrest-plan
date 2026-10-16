# 04 · 端到端流程 runbook（分步可中断执行规格）

> 目标读者：三个月后的你自己。初试结束到 3 月复试之间，每周只能挤出 2 小时。
> 本文只写**命令、文件用途、验收判据、内存/时间预算**，不贴大段源码。源码在实现阶段按本文规格写进 `code/`。
>
> 设计前提（如果这两条不成立，整份文档作废）：
> 1. **16 GB 内存是硬天花板**，而且 R 内部还有一个更硬的 16 GiB 向量堆上限（实测，见 §2.3）。
> 2. **每个 step 独立可跑、产物落盘、从断点续**。任何一步结束后 `Ctrl-C` 都不该导致重跑前面的步骤。

---

## 1. 目录结构与命名约定

在项目根 `spermatogenic-arrest-plan` 下：

```
code/                    每个 step 一个脚本，编号与 step 编号一致
  00_check_env.sh        环境自检（shell + 内嵌 R/python 探针）
  01_download.sh         下载 + 校验（调用 scripts/fetch_geo.py）
  02_prepare_expression.R
  03_de_per_dataset.R
  04_act1_convergence.R
  05_deconvolution.py    ← 主路径用 Python（流式，理由见 §3.2）
  05_deconvolution.R     ← 备用 R 路径（只在 Python venv 坏掉时用）
  06_act2_null_model.R
  07_act3_residual.R
  08_calibration.py
  09_figures_tables.R
data/
  raw/<GSE>/             原始 supp 文件，只读，下载后永不修改
  raw/_manifest.tsv      下载清单（含 SHA-256），见 §4 Step 1
  interim/<GSE>_expr.rds 每数据集单独的表达对象（List 或 dgCMatrix），**不跨数据集合并**
  interim/<GSE>_meta.tsv 样本 → 分组（WT/KO）映射，人工核对过
results/
  de/<GSE>_de.tsv        统一格式 DE 表
  de/_index.tsv          DE 表索引
  deconv/                S 特征谱、P 比例、诊断
  act1/ act2/ act3/      三幕各自的中间结果
  calib/                 校准输出
figures/                 全部图，文件名 = 图号 + 短名（fig01_...）
reports/                 本文 + 其他报告
envs/                    setup.sh（环境搭建脚本）
.Rlib/                   项目内 R 库（系统库不可写）
.venv/                   项目内 Python venv（**有 numpy/pandas/scipy/sklearn**）
```

**命名约定（强制，否则断点续跑会找错文件）**

- 数据集一律用 `GSE` 号大写，不带任何后缀。
- DE 表固定列名与顺序：`gene, logFC, AveExpr, t, P.Value, adj.P.Val`。TSV，带表头，`gene` 唯一且不重复。
- 中间产物一律 `.rds`（R）或 `.parquet`/`.npz`（Python）；**不用 `.RData`**，因为它会把整个工作区写进去，恢复时污染全局环境。
- 每个 step 结束时写一个 `<step>.done` 哨兵文件（内容 = 该步产物的 SHA-256）。续跑脚本先查哨兵，命中就跳过。

**外部工具的实际路径（本机实测）**

| 工具 | 路径 | 状态 |
|---|---|---|
| Python（分析用，有全套包） | `.venv/bin/python` | numpy 2.5.4 / pandas 3.0.6 / scipy 1.18.1 / sklearn 1.9.1 / matplotlib 3.11.2 / openpyxl 3.1.5 |
| Python（下载器用） | 任意 python3 | `scripts/fetch_geo.py` 只用标准库 |
| ⚠️ `.venvtest/bin/python` | — | **只有 numpy**，pandas/scipy 全无。别用它 |
| ⚠️ `python3`（conda base 3.14.7） | — | **没有 numpy**。别用它做分析 |
| R | `Rscript` | 4.6.1 |
| R 库 | `.libPaths(c(".Rlib", .libPaths()))` | 必须写在每个脚本第一行 |
| 下载器 | `scripts/fetch_geo.py <url> <dest>` | **6 并发 range，3 MB/块，可断点续传**（见 §2.6） |

---

## 2. 实测硬约束（这些数字都是本机量出来的，不是估的）

### 2.1 机器

| 项目 | 实测值 |
|---|---|
| CPU | Apple M4，10 核 |
| 内存 | **16.0 GB**（`hw.memsize` = 17,179,869,184 B） |
| R 向量堆硬上限 | **`mem.maxVSize()` = 16,384 Mb = 16 GiB**（这是 R 自己的墙，不是操作系统的） |
| 磁盘可用 | **82 GB**（`/System/Volumes/Data`，已用 123 GB / 228 GB） |
| conda | **不可用**，全文不出现 conda/mamba |

### 2.2 GSE112393 参考矩阵的真实尺寸（早期估算有误）

从 gzip 文件头读出的原始文件名：`mergedMouseAdultST25_34633cells_37241genes.dge.txt`

| 项目 | 值 |
|---|---|
| 基因数 G | **37,241**（早期估算 ~20k，**低估了 1.9 倍**） |
| 细胞数 N | **34,633** |
| 稠密 float64 | **10.32 GB**（= 9.61 GiB）— 早期估算 5.5 GB |
| 稠密 float32 | 5.16 GB |
| 压缩包 | 81,261,540 B = 81.3 MB |
| 解压后 | **2,584,100,928 B = 2.58 GB**（`gzip -l` 与 `gzip -dc | wc -c` 双重实测一致，已对完整文件核验） |
| 总行数 | **37,242**（1 表头 + 37,241 基因行，已实测） |
| 表头字段数 / 数据行字段数 | **34,633 / 34,634**（后者 = 基因名 + 34,633 个数值，已实测） |

**实测的失败模式**（比"内存不够"更精确）：

```r
G <- 37241; N <- 34633
m <- matrix(0, nrow=G, ncol=N)   # ✅ 成功，object.size 报 9.6 Gb
m2 <- m * 2                      # ❌ 报错：
#   cannot allocate vector of size 9.6 Gb
#   vector memory limit of 16.0 Gb reached, see mem.maxVSize()
```

含义：**这个矩阵能在 R 里"分配"出来，但一动它就死。** 任何 `t()`、`scale()`、子集复制、`as.matrix()` 都会再造一份，10.32 × 2 = 20.6 GB > 16 GiB → 直接崩。
所以"先把矩阵读进来再说"这条路**从第一步就是死的**，不是"跑慢一点"的问题。

### 2.3 各步矩阵尺寸与内存（公式：字节数 = G × N × 8）

| 对象 | G | N | float64 | float32 | 判断 |
|---|---|---|---|---|---|
| GSE112393 原始 DGE（稠密） | 37,241 | 34,633 | **10.32 GB** | 5.16 GB | ❌ 禁止 |
| 参考系基因子集 5,000 | 5,000 | 34,633 | 1.39 GB | 693 MB | ⚠️ 可行但没必要 |
| 参考系基因子集 2,000 | 2,000 | 34,633 | 554 MB | 277 MB | ✅ 安全 |
| **细胞类型特征谱 S（流式产出）** | 37,241 | 11 | **3.3 MB** | — | ✅ 主路径 |
| S 的 marker 子集 | 2,000 | 11 | 0.18 MB | — | ✅ |
| 细胞比例 P | 11 | ~100 | 9 KB | — | ✅ |
| 单个 bulk 数据集（最大 n=40） | 20,000 | 40 | 6.4 MB | 3.2 MB | ✅ |
| 单个 bulk 数据集 + DESeq2 全开销 | 20,000 | 12 | 实测峰值 **0.93 GB** | — | ✅ |

> 关键结论：**唯一危险的对象是那个 37,241 × 34,633 的稠密矩阵。** 只要不构造它，整条流水线的峰值内存都不到 2 GB。
> 而 §4 Step 5 的流式方案恰恰不需要构造它——它只需要 3.3 MB 的 S。

### 2.4 实测性能（本机，供时间预算用）

| 操作 | 规模 | 实测耗时 |
|---|---|---|
| gzip 解压吞吐 | — | 958 MB/s（解压后字节计） |
| 逐行 split + `np.asarray(float64)` 解析 | 20.9 MB | 55.8 MB/s 端到端 |
| **完整 DGE 单次流式扫描（解压+解析+按 11 类聚合）** | **2.58 GB**（37,241 基因行） | **实测 50.9 s**（50.8 MB/s） |
| 流式聚合峰值 RSS | **完整文件，含 S** | **实测 48.6 MB** |
| DESeq2 | 20,000 × 12 | 2.8 s |
| limma-voom | 20,000 × 12 | 0.6 s |
| R `nnls` | 2,000×9，80 样本 | 0.03 s |
| scipy `nnls` | 2,000×9，80 样本 | 0.02 s |

### 2.5 R 包可用性（**以本节实测为准**）

早期以为缺 edgeR/nnls/UpSetR——**实际已装好**。实测：

| 包 | 实测 | 说明 |
|---|---|---|
| limma, DESeq2, **edgeR**, clusterProfiler, org.Mm.eg.db, GEOquery, data.table, Matrix, ggplot2, pheatmap, patchwork, dplyr, readxl | ✅ OK | |
| **nnls** | ✅ OK（在 `.Rlib`） | 去卷积可用 R 路径 |
| **UpSetR** | ✅ OK（在 `.Rlib`） | |
| quadprog, sva, fgsea, GSVA, ComplexHeatmap, preprocessCore, SingleCellExperiment | ❌ MISSING | 见下 |
| **affy / illuminaio** | ❌ MISSING | 芯片 CEL/IDAT 无法从原始文件处理 |

**替代方案**：`ComplexHeatmap` → `pheatmap`（已装）；`sva` 不需要（本项目**不合并**数据集，无 batch 校正需求）；`preprocessCore`/`affy` 不需要（芯片走 `GEOquery::getGEO` 直接拿作者已归一化的 series matrix，**不要**下 CEL）；`fgsea`/`GSVA` 是可选的富集分析，不阻塞主线。

### 2.6 下载的真实情况（实测）

> **⚠️ 本节已按 `scripts/fetch_geo.py` 的当前版本（166 行）重写。** 早期版本的分块先攒在内存 `parts` 字典里、中断即全丢；那个版本**已被替换**。下面写的是**现在磁盘上这个版本**的行为，已实测验证。

**现在的 `fetch_geo.py` 是可断点续传的**（实测）：
- `CHUNK = 3_000_000`（3 MB）、`STREAMS = 6`、`RETRIES = 6`、`max_rounds = 8`。
- `scripts/fetch_geo.py:74` `partdir = dest + ".parts"`，块文件命名 `part_{i:05d}`。
- `scripts/fetch_geo.py:105-108` **每块到达后立刻写自己的 `.tmp` 再 `os.replace` 成正式块** —— 不再攒在内存里。
- `scripts/fetch_geo.py:83-99` 每轮开始前扫描已有块，`os.path.getsize(p) == want` 的跳过 → **重跑即续传**。
- `scripts/fetch_geo.py:132-137` 组装前校验每块尺寸，缺块则报 `still missing N chunk(s)` 并提示重跑，**不会产出半截的成品文件**。
- 内存上界 = `streams × chunk` = 6 × 3 MB = **18 MB**，与文件总大小无关（旧的"GB 级文件打爆内存"问题已不存在）。

**实测验证（4 块 × 400 KB 的小文件，强制 `--chunk 400000`）**
| 测试 | 结果 |
|---|---|
| 完整下载后与原始文件 `cmp` | ✅ 字节完全一致 |
| 删掉 `part_00002` 后重跑 | 只重下 1 块（`1 chunk(s) to fetch (1.1/1.5 MB already on disk)`），成品仍字节一致 |
| 把 `part_00001` 截成 0 字节后重跑 | 被正确识别为缺失并重下，成品一致（**尺寸校验有效**） |
| HEAD 404 的文件 | 报 `cannot determine size` 且**不创建目标文件**（不产生垃圾文件） |
| 实测吞吐 | 1.5 MB / 约 20-29 s ≈ **50-72 KB/s**（小块、含 503 退避） |

**坑 1 —— NCBI 会返回 HTTP 503 限流，且会成片持续。** 实测 8 并发列目录时大面积 503；改成串行 + 指数退避后成功。
→ **对策：列目录/探尺寸串行 + 重试；只有下正文才开 6 并发。**

**坑 2 —— 503 返回的是 HTML 错误页，而 curl/urllib 不会把它当错误。** 写本文时复测：对同一个 supp 文件连续 6 次请求（退避 15/30/45/60/75/90 s）**全部返回 503 HTML**，`gzip -dc` 报 `not in gzip format`；而与此同时对另一个文件的单次请求却成功。特征：
- 正文是 `<!DOCTYPE html>...<title>Service unavailable!</title>`，状态码 503，**`curl -s` 下静默写出 HTML 文件且 exit code 仍为 0**。→ **必须检查 magic bytes，不能只看 exit code 或"文件非空"。**
- 同一时间不同文件的命运不同 → 是按路径/连接限流，不是全站宕机。
- 退避到 90 s 仍未通过 → **长退避未必能救**。中断 10-20 分钟后再试，或换 `https://www.ncbi.nlm.nih.gov/geo/download/?acc=...&format=file&file=...` 这条路径，往往立刻就好。
- **好消息**：当前版 `fetch_geo.py` 用 `urllib` 且校验 `len(data) == want`（`scripts/fetch_geo.py:60-62`），503 的短正文会被判为 `short read` 并进入重试，不会污染 `part_*`。所以**经过 `fetch_geo.py` 下载的文件天然免疫这个坑**；但**用 `curl`/`wget` 手工下载时必须自己查 magic bytes**。

**必备的下载校验（每个文件都做，别信 exit code）**

```bash
# 1) 确认是 gzip 而不是 503 HTML（前 2 字节必须是 1f 8b）
head -c 2 "$f" | xxd -p          # 期望输出 1f8b；若为 3c21（<!）就是 503 页
# 2) 确认不是 503 正文
head -c 200 "$f" | grep -qi 'service unavailable\|<!DOCTYPE' && echo "503! 删掉重下"
# 3) 尺寸必须等于 HEAD 的 Content-Length
[ "$(stat -f%z "$f")" = "$expected" ] || echo "尺寸不符"
# 4) 只有过了 1-3 才写进 _manifest.tsv
```

> 这条**否定了"重试就能成功"的乐观假设**：Step 1 必须把"确认是 gzip magic"写成硬性验收，否则你会在 Step 5 才发现拿到的是一份 HTML 错误页，白白浪费一个会话。

> **续传语义的两个边界（实现 Step 1 时必须知道）**
> 1. **续传的粒度是块，不是字节。** 一个 3 MB 的块若只下到 2.9 MB，重跑会**整块重下**（尺寸不符即丢弃重来），不会从 2.9 MB 接着下。所以最坏浪费 3 MB，可接受。
> 2. `.parts` 目录**不会自动清理**，成功组装后仍保留（`scripts/fetch_geo.py:150-151` 提示可 `rm -rf`）。参考系的 `.parts` 是 27 × 3 MB ≈ 78 MB，**收工前记得删**，否则会悄悄吃掉磁盘。

**实测文件尺寸表**（`HEAD` 量的，用来算下载预算）：

| GSE | 文件 | 大小 | 是否要下 |
|---|---|---|---|
| GSE112393 | `MergedAdultMouseST25_DGE.txt.gz` | 81.3 MB | ✅ **必需**（参考系） |
| GSE112393 | `MergedAdultMouseST25_PerCellAttributes.txt.gz` | 1.5 MB | ✅ **必需**（细胞注释） |
| GSE112393 | `RAW.tar` | 107.8 MB | ❌ 不要（与 DGE 重复） |
| GSE112393 | `12GermCellClusters_AllGeneExp.xlsx` | 11.4 MB | ❌ 不要（粒度不对） |
| GSE236388 | `XAP5L_counts.txt.gz` + `XAP5_counts.txt.gz` | 0.8 MB | ✅ 计数矩阵 |
| GSE213442 | `cpm_matrix.txt.gz` | 1.4 MB | ✅ |
| GSE214201 | `cpm_matrix.txt.gz` | 1.0 MB | ✅ |
| GSE318472 | `FPKM_Spry4.csv.gz` | 1.0 MB | ✅（FPKM，非计数） |
| GSE221433 | `RAW.tar` | 7.9 MB | ⚠️ 需解包 |
| GSE235777 | `RAW.tar` | 2.0 MB | ⚠️ 需解包 |
| GSE259417 | `P16_DEG.csv.gz` | 2.4 MB | ⚠️ **只是作者算好的 DEG，不是计数矩阵** → 能进第一幕（有 logFC），**不能进去卷积** |
| GSE70889 | 2 个 `read_group_tracking*.txt.gz` | 4.6 MB | ⚠️ 分选肾母细胞 |
| GSE55180 | `RAW.tar` | 11.7 MB | ⚠️ 需解包 |
| GSE121326 | `RAW.tar` | 8.5 MB | ⚠️ 芯片 CEL |
| GSE158471 | `RAW.tar` | 1.1 MB | ⚠️ 需解包 |
| GSE294209 | `RAW.tar` + geneTPM + isoTPM | 16.2 MB | ✅ 有 TPM |
| GSE137307 | `RAW.tar` | **540.1 MB** | ❌ 芯片 CEL，且 affy 未装 |
| GSE253154 | `RAW.tar` + DEG xlsx | **2,670.9 MB** | ❌ 除非必要时再下 |
| **GSE244681** | `RAW.tar` | **496,540 MB（496 GB）** | ❌❌❌ **绝对不要下**，磁盘只有 82 GB，而且下载要 64,000 分钟 |

> 参考系下载预算：81.3 + 1.5 = **82.8 MB**，按 126 KB/s = **约 11 分钟**（单流要 70+ 分钟，所以必须用 `fetch_geo.py`）。

**入选数据集清单以 `data/dataset_selection.csv` 为准。** 本文 Step 1 的清单模板从那个文件生成；如果它还不存在，用上面的表 + `data/inventory.csv` 里 `usability=可直接用 且 tissue=全睾丸/待查` 的行生成初版。

---

## 3. 内存控制总策略（先定原则，后面每步引用）

### 3.1 铁律

1. **永远不构造 37,241 × 34,633 的稠密矩阵**。这是唯一会撑死的对象，任何形式的"先读进来"都禁止。
2. **矩阵字节数预算**：任何时刻所有大对象之和 ≤ 6 GB（留 10 GB 给 R 的临时副本和 OS）。
3. 大对象用完立刻 `rm(x); gc()`（R）或 `del x; gc.collect()`（Python）。
4. 跨 step 传数据一律落盘，不做内存里的管道。

### 3.2 参考系为什么用 Python 流式，而不是 R

- R 侧 `data.table::fread` 读一个 3.75 GB 的 TSV 会先整体物化，然后还要 `as.matrix()` → 两份 10.32 GB → 必死。
- Python 侧可以**逐行处理，永不物化**：每读一个基因行（37,241 行之一），得到一个长度 34,633 的向量（277 KB），除以该细胞的 size factor，按细胞类型 `bincount` 累加进 S 的一行，然后丢弃这个向量。峰值 = 单行缓冲 + S（3.3 MB）+ 全基因组 S_final。
- 实测峰值 RSS **47.7 MB**。

### 3.3 为什么不需要"分块读大矩阵"这种复杂方案

直觉上会想"那就分块读 2000 个基因的块"。但实测表明有更简单的解法：**伪 bulk 聚合不需要细胞 × 基因矩阵，只需要聚合结果 S（基因 × 细胞类型）**。
对比：

| 方案 | 峰值内存 | 需要几趟扫描 | 复杂度 |
|---|---|---|---|
| 稠密读入 | 10.32 GB ×2 | 1 | 崩 |
| 基因分块（每块 1000 行，f4） | 139 MB/块 + 拼接 | 1 | 中 |
| **流式聚合到 S** | **< 300 MB** | **2**（见 §4 Step 5） | 低 |

两趟扫描 = 2 × 67 s = **约 2.3 分钟**，代价可忽略。**采用流式聚合。**

---

## 4. 分步规格

每一步统一给：**输入 / 输出 / 命令 / 验收判据 / 耗时 / 峰值内存 / 可中断点 / 失败怎么办**。

---

### Step 0 · 环境自检

| 项 | 值 |
|---|---|
| 脚本 | `code/00_check_env.sh`（把 R 和 Python 探针都塞进一个 shell，避免两个文件） |
| 输入 | 无 |
| 输出 | `reports/00-env-report.txt`（体检表）+ 退出码（0 = 全绿） |
| 耗时 | **< 30 s** |
| 峰值内存 | < 500 MB |
| **可中断点** | 无意义（30 秒的事，跑完为止） |

**检查项（逐条打印 PASS/FAIL）**

1. 工作目录正确（`pwd` 末尾是 `生精阻滞共性分析`）。
2. `.Rlib` 存在；`Rscript -e '.libPaths(c(".Rlib",.libPaths())); rownames(installed.packages())'` 能列出包。
3. **关键包逐个探测**：`limma, DESeq2, edgeR, nnls, GEOquery, data.table, Matrix, ggplot2, pheatmap, patchwork, dplyr, org.Mm.eg.db, UpSetR` → 期望全 OK。缺失的打印名字，不要中止。
4. `.venv/bin/python` 能 `import numpy, pandas, scipy, sklearn, matplotlib` 并打印版本。
5. **显式警告 `.venvtest` 和 conda base python 不可用于分析**（它们缺 numpy/scipy）。
6. `mem.maxVSize()` 打印出来，并断言 `== 16384 Mb`（如果 R 版本变了，这个数会变，需要重新评估 §3）。
7. 内存：`sysctl -n hw.memsize` → 期望 16 GB。
8. 磁盘：`df -h` 至少 60 GB 可用。
9. 网络：对 `https://ftp.ncbi.nlm.nih.gov/` 做一次 `HEAD`，**带 3 次重试和退避**（实测会 503）；报告耗时和等效吞吐。
10. `scripts/fetch_geo.py` 存在且 `python3 -c "import ast;ast.parse(open('scripts/fetch_geo.py').read())"` 通过。
11. **断言禁忌项**：`which conda mamba` 如果命中，打印警告说本项目不使用 conda（防三个月后自己手滑）。

**验收判据**：体检表里 `.Rlib` 关键包全 OK、`.venv` 五个包全 OK、磁盘 ≥ 60 GB、`mem.maxVSize == 16384`。网络项允许 FAIL（离线也能干活，Step 1 单独处理）。

**如果这步失败怎么办**
- `.Rlib` 关键包缺失 → 回去跑 `envs/setup.sh`。**不要** `install.packages` 到默认库（不可写）。
- `.venv` 坏 → `python3 -m venv .venv && .venv/bin/pip install numpy pandas scipy scikit-learn matplotlib openpyxl`。注意用 `python3 -m venv`，**不要** conda。
- `mem.maxVSize` 不是 16384 → 说明 R 版本变了；重跑一遍 §2.2 那个"分配 + 复制"的探针，确认新的上限，回来改 §3 的 6 GB 预算。
- 网络失败 → 不阻塞，直接进 Step 1，Step 1 自己会重试。

---

### Step 1 · 下载与核验

| 项 | 值 |
|---|---|
| 脚本 | `code/01_download.sh`（读清单 → 循环调用 `scripts/fetch_geo.py` → 写 manifest） |
| 输入 | `data/dataset_selection.csv`或 `data/inventory.csv` |
| 输出 | `data/raw/<GSE>/...`、`data/raw/_manifest.tsv` |
| 耗时 | 参考系 **11 min**；加 bulk 数据集按 §2.6 表累加（多数 < 2 min/个） |
| 峰值内存 | 单文件分块缓存：`ceil(size/2MB) × 2MB`。81 MB 文件 → 82 MB。**> 300 MB 的文件必须先确认**（§2.6 坑 1） |
| **可中断点** | ★ **每个文件下完后**。`_manifest.tsv` 每次追加写入并 `fsync`。续跑时跳过已校验通过的文件 |

**下载清单模板** `data/download_plan.tsv`（由 Step 0 之后手工生成）：

```
gse	url	size_bytes	sha256	needed_for	priority
GSE112393	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112393/suppl/GSE112393_MergedAdultMouseST25_DGE.txt.gz	81261540	<下载后填>	deconvolution_reference	1
GSE112393	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112393/suppl/GSE112393_MergedAdultMouseST25_PerCellAttributes.txt.gz	1495253	<下载后填>	deconvolution_reference	1
GSE236388	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236388/suppl/GSE236388_XAP5L_counts.txt.gz	<HEAD>	<填>	bulk_de	2
```

**命令形态**（逐个文件，串行启动以避免 503）：

```bash
# 1) 先串行探尺寸并写清单（带退避重试，应对 503）
#    用途：确认每个文件多大，避免误下 496 GB 的 GSE244681_RAW.tar
python3 code/01_download.sh --probe-only

# 2) 人工过一眼清单，删掉体积不合理的行（尤其 RAW.tar）
less data/download_plan.tsv      # 用途：确认没有 >300 MB 的非必要文件

# 3) 逐个下载。每个文件 6 并发 range，3 MB/块，可续传。
#    中断了就直接重跑同一条命令——已完成的分块会被跳过（实测）。
while IFS=$'\t' read -r gse url size sha rest; do
  dest="data/raw/$gse/$(basename "$url")"
  [ -f "$dest" ] && continue                                  # 成品已存在，跳过
  mkdir -p "data/raw/$gse"
  python3 scripts/fetch_geo.py "$url" "$dest"                 # 可反复重跑续传
done < data/download_plan.tsv

# 4) 校验：gzip magic + 尺寸 + SHA-256，写 manifest
for f in data/raw/*/*/*.gz; do
  [ "$(head -c 2 "$f" | xxd -p)" = "1f8b" ] || echo "❌ 非 gzip（503 页？）: $f"
done
shasum -a 256 data/raw/*/*/*.gz > data/raw/_sha256.txt        # 用途：留档，后续引用
python3 code/01_download.sh --verify                          # 用途：对比清单尺寸，标出不匹配

# 5) 收工清理分块目录（不删会一直占磁盘）
find data/raw -maxdepth 2 -name '*.parts' -type d -exec rm -rf {} +
```

**验收判据**
- 每个文件的 `stat -f%z` == 清单里的 `size_bytes`（`fetch_geo.py` 自己会查一次，这里查第二次）。
- `data/raw/_manifest.tsv` 行数 == 下载清单行数，且每行都有 64 位 hex 的 SHA-256。
- 参考系的 SHA-256 与 数据集报告里记录的一致。
- ~~特别注意 `GSE112393_DGE.txt.gz` 当前是 0 字节~~ → **已过时（2026-10-16 更正）**：该文件现已完整下载，`data/raw/ref/GSE112393_DGE.txt.gz` = **81,261,540 B**，`gzip -t` 通过。`fetch_geo.py` 已加固：下载前先 `HEAD` 比对尺寸 + `gzip_ok()` 全量校验 CRC/ISIZE，已完整则打印 `already complete: ... -- skipping` 直接跳过。**参考系不再是阻塞点。**

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| HTTP 503 | 退避重试。列目录/探尺寸串行；正文交给 `fetch_geo.py`（6 并发 + 每块 6 次退避重试，503 短正文会被判 `short read` 自动重来） |
| 中断后 `.part` 残留 | `fetch_geo.py` **已支持块级续传**（见 §2.6）：重跑即跳过已完整的块。**不要手动删 `.parts`** —— 拼装成功后脚本自己会 `rmtree`。若某块半途断在 `.tmp`，脚本启动时会自动清理。 |
| 文件 > 300 MB | **先停下来问自己是否真的需要**。RAW.tar 里往往混着 ChIP/Hi-C，解包后 90% 用不上 |
| SHA 不匹配 | 删文件重下。连续两次不匹配 → 检查是不是被 GEO 重新发布过（用 `HEAD` 看 `Last-Modified`） |
| 参考系 DGE 下不动 | **已不再阻塞**：该文件已完整下载（81,261,540 B，`gzip -t` 通过）。若将来重下遇到持续 503，降级：先用 `PerCellAttributes`（1.5 MB）做注释核对，DGE 挂后台续传 |

---

### Step 2 · 每数据集单独读入 + 统一基因空间

| 项 | 值 |
|---|---|
| 脚本 | `code/02_prepare_expression.R` |
| 输入 | `data/raw/<GSE>/*`、`data/raw/_manifest.tsv` |
| 输出 | `data/interim/<GSE>_expr.rds`、`data/interim/<GSE>_meta.tsv`、`data/interim/02.done` |
| 耗时 | 每数据集 **5-30 s**（读 + 注释）；8 个 ≈ 2-4 min |
| 峰值内存 | **< 1.5 GB**（20,000 × 40 f8 = 6.4 MB，主要开销是解包和 GEOquery） |
| **可中断点** | ★ **每个数据集处理完后**。`02.done` 里记已完成 GSE 列表 |

**做四件事**

1. **读入**：计数矩阵（`GSE236388_*_counts.txt.gz`、`GSE213442_cpm_matrix.txt.gz` 等）用 `data.table::fread`；`RAW.tar` 先 `tar -tf` 列出内容再只解需要的成员（**不要整包解**）；芯片走 `GEOquery::getGEO` 拿系列矩阵（**不要下 CEL，affy 没装**）。
2. **统一基因符号**（每数据集**内部**做，绝不跨数据集对齐数值）：
   - `org.Mm.eg.db` 做 `SYMBOL` 规范化（大小写统一成首字母大写、其余小写；Ensembl ID → symbol；旧别名 → 现名）。
   - 一个 symbol 映射到多个行（多 Ensembl）→ 按行和合并，并记录合并日志到 `data/interim/<GSE>_gene_collapse.tsv`。
   - 映射不上的 → 丢弃并计数，写进日志（正常丢 5-15%）。
   - 版本后缀（`Gm12345.1`）剥掉。
3. **过滤**（在每个数据集**内部**用统一阈值，但不同数据集基因数可以不同）：
   - 低表达：`rowSums(counts >= 10) >= min(3, n_min_group)`（至少 3 个样本或该数据集最小分组样本数）——只有 2 样本/组的（如 GSE108993 n=5 分成 3/2）放宽到 ≥ 2。
   - 线粒体：`^mt-` 开头的全部剔除，并记录剔除前后的总表达占比（这是组成变化的**诊断量**，Step 6 会用到）。
   - 核糖体、性染色体基因**不要**剔（`%ChrX`/`%ChrY` 在参考系注释里就有，是重要协变量）。
4. **样本元数据人工核对**：`data/interim/<GSE>_meta.tsv` 列 `sample, group(WT|KO), node, tissue, platform, n_replicates`。`group` 必须手工看过 `sample_titles` 才准（`data/inventory.csv` 的 `n_ko`/`n_wt` 有已知不一致，例如 GSE235777 记 3/6）。

**绝对禁止**：把不同数据集的表达值拼成一个矩阵，或对拼起来的矩阵做 `removeBatchEffect`/`ComBat`。本项目**只在 logFC 层面比较**。在脚本里加一句断言：任何 `cbind`/`merge` 跨 GSE 的操作直接 `stop()`。

**验收判据**
- 每个 `_expr.rds` 里 `dim()` 报的基因数 ≥ 12,000（低于这个数说明基因符号映射出了问题）。
- `_meta.tsv` 的 `group` 只有 `WT`/`KO` 两个取值，且两组的样本名与 `data/inventory.csv` 的 `sample_titles` 对得上。
- 每个数据集 `colSums` 无 0（有 0 说明读错了朝向——**矩阵转置错误是最常见的读入 bug**，务必打印 `dim()` 和 `head(rownames())` 人工看一眼）。
- `02.done` 里 8 个 GSE 全部出现。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| 基因数 < 12,000 | 多半是 symbol 映射成了 Ensembl 或大小写没统一。打印未映射的前 20 个 ID 看看 |
| 矩阵朝向搞反 | 用 `dim()` 判断：行数应该 ≈ 20,000（基因），列数 ≈ 样本数。反了就 `t()` |
| `RAW.tar` 里是 `BW`/`bigwig` 而非计数 | 该数据集**放弃**（如 GSE227944 是 `signal_only`），从清单里划掉，不要硬凑 |
| GSE259417 只有作者 DEG | 跳过 Step 2/3，直接从它的 DEG 提取 logFC，标 `de_source=author`，**并标记它不能进去卷积（Step 5 跳过）** |
| FPKM 而非计数（GSE318472） | 可以进 limma（`lmFit` 直接吃 log2(FPKM+1)），但**不能进 DESeq2**。在元数据里标 `quant=FPKM` |
| 内存不够 | 几乎不可能。真发生了说明 `fread` 读了畸形文件，加 `--probe` 看前 5 行 |

---

### Step 3 · 每数据集单独 DE

| 项 | 值 |
|---|---|
| 脚本 | `code/03_de_per_dataset.R` |
| 输入 | `data/interim/<GSE>_expr.rds` + `_meta.tsv` |
| 输出 | `results/de/<GSE>_de.tsv`（**固定列** `gene, logFC, AveExpr, t, P.Value, adj.P.Val`）、`results/de/_index.tsv`、`results/de/03.done` |
| 耗时 | 计数 + limma-voom **0.6 s** / DESeq2 **2.8 s**（实测 20,000×12）；8 个数据集含 I/O ≈ **1-2 min** |
| 峰值内存 | **< 1.2 GB**（实测 DESeq2 单数据集峰值 0.93 GB） |
| **可中断点** | ★ **每个数据集跑完后**。`_index.tsv` 逐行追加 |

**方法选择（按 `quant` 字段自动分派）**

| 数据 | 方法 | 理由 |
|---|---|---|
| 原始计数（counts） | `DESeq2`（主）+ `limma-voom`（敏感性） | DESeq2 对小 n 更稳；两个都跑能验证结论不依赖方法 |
| 原始计数且 n ≥ 8 | `limma-voom`（快，0.6 s） | 大 n 时两者结果高度一致 |
| FPKM/TPM/CPM | `limma::lmFit` on `log2(x+1)` | 非计数不能进 DESeq2 |
| 芯片 | `limma::lmFit` on 系列矩阵 | 已归一化，直接线性建模 |
| 作者已算 DEG | 直接解析，重命名列 | 标 `de_source=author` |

**设计要求**
- 设计矩阵：`~ group`（两组）。有配对的（同窝、同批）加 `~ batch + group`，但**只有当元数据里真有 batch 字段才加**。
- 单样本/组的极端情况（GSE108993 是 3/2、GSE174731 各 3 只混池算 2 样本）→ DESeq2 的 dispersion 估计会不稳，**强制加 `apeglm` 或 `ashr` 收缩，并在 `_index.tsv` 标 `unstable_ok=FALSE`**，Step 4 里做敏感性剔除。
- 全部基因输出（`number=Inf, sort.by="none"`），不要只输出显著基因——第二幕需要全基因组的 logFC 向量。
- `p.adjust` 用 BH。芯片用 `limma::topTable` 自带的 `adj.P.Val`。
- 每个数据集额外写一行诊断到 `_index.tsv`：`n_genes, n_sig(padj<0.05), n_sig(|logFC|>1), frac_mt_before, frac_mt_after, method, de_source, unstable`。

**验收判据**
- 8 个 `_de.tsv` 全部存在，列名与顺序完全一致（脚本末尾跑一个 schema 断言）。
- 每个表的 `gene` 列无重复，`logFC` 无 `NA`（`NA` 说明设计矩阵奇异，必须查）。
- 抽查：敲除模型里 `Prm1`/`Tnp1` 应显著下调，`Sox9`/`Wt1` 应上调。**如果没看到这个方向，立刻停下**——说明分组标签反了（WT/KO 搞反是第二常见的 bug）。
- `_index.tsv` 行数 == 数据集数。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| `logFC` 全 NA | 设计矩阵奇异，多半是某组只有 1 个样本或分组全一样。检查 `table(group)` |
| Prm1 方向反了 | WT/KO 标签反了。回 Step 2 改 `_meta.tsv` |
| DESeq2 报 dispersion 错误 | 降到 limma-voom（`edgeR` 已装） |
| 老芯片数据读进来是 log 值 | `limma` 直接吃，不要 `log2` 两次 |
| 某个数据集 n=2 且分散 | 保留但标 `unstable`，Step 4 做剔除敏感性 |

---

### Step 4 · 第一幕：跨模型 logFC 相似度（天真收敛）

| 项 | 值 |
|---|---|
| 脚本 | `code/04_act1_convergence.R` |
| 输入 | `results/de/*_de.tsv` |
| 输出 | `results/act1/` 下的相似度矩阵（`.tsv`）、聚类树（`.rds`）、`act1/04.done` |
| 耗时 | **< 1 min**（8 个 20,000 长的向量，相关矩阵是 8×8） |
| 峰值内存 | **< 500 MB** |
| **可中断点** | ★ 跑完即止（1 分钟），无需中途断点 |

**做法**
1. **取共同基因**：所有数据集 logFC 向量取交集（用统一后的 symbol。跨平台的话交集会缩到 ~10,000-13,000，这是预期内的，记录下来）。
2. **两套相似度**，都要：
   - **Spearman 秩相关**（主）：跨平台可比，对平台差异不敏感。
   - **Pearson on logFC**（辅，仅 RNA-seq 子集）：数值可比时用。
   - **Top-N 重叠**（Jaccard / 超几何）：对 top 100/500/2000 各算一次，看结论稳不稳。
3. **聚类**：对 `1 - Spearman` 距离矩阵做 `hclust(method="average")`，`cutree` 到 3 组，把聚类结果与 `node` 标注对照（这是项目笔记 5.1 说的"可验证的预测"）。
4. **敏感性**：把 `unstable=TRUE` 和 `de_source=author` 的数据集各剔除一次，重跑相似度，看矩阵是否稳定。三个版本都存。

**验收判据**
- 相似度矩阵对角全 1，对称。
- 平均 off-diagonal Spearman **预期 > 0.5**（这是"复现文献说法"的关键；如果只有 0.2，说明第一幕的前提不成立，整个故事要重新想——这本身是个重要发现，必须如实记录）。
- 聚类结果与 `node` 的一致性能用 ARI 或简单的对照表量化并写进 `results/act1/node_vs_cluster.tsv`。
- **图先不画**，Step 9 统一画（避免重复劳动）。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| 共同基因 < 5,000 | 多半是 symbol 映射不统一。回 Step 2 检查映射规则是否每个数据集都一致 |
| 平均相关 < 0.2 | **不要掩盖**。这直接推翻第一幕的前提。记录之，并与零模型规格讨论是否需要重构故事 |
| 某个数据集把整体拉偏 | 看相似度矩阵里它的那一行，单独剔除后重算，作为敏感性版本 |

---

### Step 5 · 去卷积（**内存风险最高的一步**）

| 项 | 值 |
|---|---|
| 脚本 | `code/05_deconvolution.py`（主路径）；`code/05_deconvolution.R`（备用） |
| 输入 | `data/raw/GSE112393/GSE112393_MergedAdultMouseST25_DGE.txt.gz`、`..._PerCellAttributes.txt.gz`、各 bulk 表达矩阵（Step 2） |
| 输出 | `results/deconv/S_celltype.tsv`（37,241×11，3.3 MB）、`results/deconv/P_proportions.tsv`、`results/deconv/markers.tsv`、诊断图、`deconv/05.done` |
| 耗时 | 2 趟流式扫描 × 67 s = **2.3 min** + marker 选择 ~10 s + NNLS ~1 min → **约 5 min** |
| 峰值内存 | **实测 47.7 MB**（流式原型）。含 pandas 开销保守估 **< 500 MB** |
| **可中断点** | ★★ **三处**：(a) 第 1 趟列和算完后；(b) 第 2 趟 S 产出后；(c) NNLS 完每个 bulk 数据集后 |

**参考系事实（实测，不是假设）**
- 文件头行是 34,633 个细胞条码（形如 `ST1_ATCTTGCACATC`），之后每行是「基因名 + 34,634 个 tab 分隔字段」。
- `PerCellAttributes` 的**精确结构（实测字节级核对，早期版本说的"第 3 行起"和本文早期版本的"第 5 行是表头"都不对）**：
  - 第 1 行 = `###### Per-cell Attributes for All 34,633 Cells...`（注释，1 列）
  - 第 2 行 = `### Cell type: 1-InnateLymph; ...`（注释，1 列）
  - 第 3 行 = `### Each row is one cell`（注释，1 列）
  - **第 4 行 = 表头**，以 `#CellBarcode` 开头，27 列
  - **第 5 行起 = 数据**，共 **34,633 行**，每行 27 列
  - 第 1 列 `#CellBarcode`、第 2 列 `Dataset`、第 3 列 `CellType`（数字 1-11）
- 注释图例（文件第 2 行给出）：`1-InnateLymph 2-Macrophage 3-Endothelial 4-Myoid 5-Leydig 6-Sertoli 7-Unknown 8-SPG 9-Scytes 10-STids 11-Elongating`
- 实测细胞类型分布（合计 34,633）：`1:64, 2:139, 3:179, 4:49, 5:314, 6:2131, 7:2205, 8:4239, 9:8792, 10:9923, 11:6598`
- 实测 Dataset 分布含 25 个来源（合计 34,633）→ **批次字段，必须当协变量**：
  `ST1:2324 ST2:1642 ST3:2720 ST4:2363 ST5:1601 ST6:2341 ST7:1549 ST8:1286 SPG1:957 SPG2:937 SPG3:1212 INT1:819 INT2:237 INT3:1210 INT4:653 INT5:1480 INT6:1453 SER1:81 SER2:161 SER3:451 SER4:716 SER5:2015 SER6:1399 SER7:3047 SER8:1979`

**⚠️ 三个实测的解析陷阱（每个都会静默产生错数据）**

1. **表头本身以 `#` 开头。** 所以"跳过所有 `#` 开头的行、把第一个非 `#` 行当表头"这个写法会**把第一条细胞 `ST1_ATCTTGCACATC` 当成表头**，然后少一个细胞、列名全是数据。正确写法是**显式跳过前 3 行**（或匹配 `^#CellBarcode` 定位表头），不要用"第一个非注释行"的启发式。
2. **文件是 CRLF 行尾**（实测 34,636 个 `\r\n`）。用 `open(...).read().split('\n')` 再 `split('\t')` 时，**每行最后一个字段会带 `\r`**（第 3 行实测为 `### Each row is one cell\r`）。读第 3 列 `CellType` 不受影响，但如果后面要读 `%MT` 之类的尾部数值列，必须 `strip()`。
3. **文件末尾没有换行符**（实测末字节为 `NA`，不是 `\n`）。后果：`wc -l` 报 **34636**、`awk 'NR>5'` 报 **34632**，**都比真实的 34,633 少 1**。如果你按 `wc -l` 或 awk 计数做断言，会得到一个幻影 off-by-one，并可能误判"少了一个细胞"。**唯一可信的计数是 `len(labels)` 或 `sum(CellType counts)`**。同理，尾部若干列含字面量 `NA`（实测末行末尾 6 个字段全为 `NA`），pandas 读的时候要显式声明 `na_values=['NA']`，否则会被当成字符串列。

**内存控制策略（这是本 runbook 最要紧的一节）**

🚫 **禁止**：`pandas.read_csv(...).values`、`np.loadtxt`、`R: as.matrix(fread(...))`、任何"把 DGE 完整读成矩阵"的写法。实测这会产生 10.32 GB 对象，R 里连一次复制都做不到（`vector memory limit of 16.0 Gb reached`）。

✅ **采用：流式伪 bulk 聚合（一趟扫描即可，永不物化）**

> **下面这段是实测跑通的版本**（在完整 2.58 GB 文件上真跑过：50.9 s、峰值 RSS 48.6 MB、S 形状 (37241, 11) 正确）。**其中有一个会直接崩溃的坑，已在下文标红**——我按最初写法跑的时候真的崩了。

```python
import gzip, numpy as np

N = 34633
# --- 1) 读细胞类型标签（注意：跳过前 4 行，不是前 3 行）---
labels = []
with gzip.open(".../GSE112393_PerCellAttributes.txt.gz", 'rt', newline='') as fh:
    for i, line in enumerate(fh):
        if i < 4: continue                       # 第 1-3 行注释 + 第 4 行表头
        f = line.rstrip('\r\n').split('\t')      # 必须 rstrip('\r')：该文件是 CRLF
        if len(f) < 3: continue
        labels.append(int(f[2]))                 # 第 3 列 = CellType
labels = np.asarray(labels, dtype=np.int64)      # 长度断言 == 34633
assert labels.size == 34633, labels.size

# ⚠️⚠️ 崩溃点：标签取值是 1..11，np.bincount 会把最大值当上界，
#      返回长度 12（index 0..11）的数组。直接塞进 (11,) 的列会报：
#      ValueError: could not broadcast input array from shape (12,) into shape (11,)
counts = np.bincount(labels, minlength=labels.max() + 1)   # -> 长度 12
assert counts.size == 12 and counts[0] == 0                # index 0 恒为 0，必须丢弃

# --- 2) 一趟扫完：同时累积 col_sum 与按类型的和 ---
col_sum = np.zeros(N, dtype=np.float64)          # 277 KB
S = np.zeros((37241, 12), dtype=np.float32)      # 12 列，第 0 列是占位、最后丢弃
ngene = 0
with gzip.open(".../GSE112393_MergedAdultMouseST25_DGE.txt.gz", 'rt') as fh:
    fh.readline()                                # 第 1 行 = 34,633 个细胞条码
    for line in fh:
        parts = line.rstrip('\n').split('\t')
        if len(parts) != N + 1:                  # 基因名 + 34,633 个数值
            continue
        x = np.asarray(parts[1:], dtype=np.float64)   # 34,633 × 8 = 277 KB 临时
        col_sum += x
        S[ngene, :] = np.bincount(labels, weights=x, minlength=12)
        ngene += 1
assert ngene == 37241, ngene

S = S[:, 1:]                                     # 丢掉 index 0 → (37241, 11)
S /= counts[1:]                                  # 除以每类细胞数 → 均值谱
np.save("results/deconv/_colsum.npy", col_sum)   # 277 KB，落盘 → 可中断点
```

**两个必须知道的实现细节**
1. **`np.bincount` 的 `minlength` 必须是 12，不是 11。** 标签是 **1..11**（不是 0..10），`bincount` 的结果长度 = `max(label)+1`。写 `minlength=11` 会得到长度 12 的数组，赋给 `S[i, :]`（形状 11）时抛 `could not broadcast input array from shape (12,) into shape (11,)`。**这是我在真文件上踩到的坑，不是推测。** 正确做法：分配 12 列、用完丢掉第 0 列。
2. **不要用 `np.fromstring(row, sep='\t')`**（早期草稿的写法）。它在 NumPy 2.x 中已 deprecated，且容易对分隔符处理不一致。用 `np.asarray(parts[1:], dtype=np.float64)`（先 `split` 再转换）——实测 50.8 MB/s；`np.fromstring` 只快约 10%，不值得冒版本风险。

**为什么一趟就够**（原稿说"两趟"是多余的）：上面这版在**同一趟**里既累加了 `col_sum` 又累加了按类型的和。归一化放在循环之后（`S /= counts[1:]`），因为按类型求和的**均值**天然就是"每类细胞的平均表达"，不需要先知道每个细胞的 size factor。**实测一趟 50.9 s 完成全部工作。**
若要按细胞做 UMI 归一化再聚合（另一种口径），才需要两趟（第 1 趟求 `col_sum`，第 2 趟除以 size factor 后 bincount），代价是 2 × 51 s ≈ 1.7 min。**两种口径都建议跑，作为敏感性分析。**

总峰值 = 单行缓冲 277 KB + S 1.8 MB + labels 277 KB ≈ **48.6 MB（完整文件实测）**

**批次效应处理**：25 个平台/批次。做两版 S：
- `S_global`：所有细胞一起聚合（主）。
- `S_within_dataset`：按 `Dataset` 分别聚合，然后取各数据集的**中位数谱**做稳健合成（敏感性）。
两版都存，Step 6 对比。

**marker 基因选择**（去卷积只需要 marker，不需要全部 37,241 个）
- 从 `S` 里按细胞类型特异性选：`特异性 = S_g,k / (mean over k of S_g,k)`，每类型取 top 200，合并去重 → 约 1,500-2,500 个 marker。
- 与 Step 2 各 bulk 数据集的基因取交集。
- **写 `markers.tsv`**（gene, cell_type, specificity），供 Step 6/7 复用。

**细胞类型粒度**：11 类里 `7-Unknown`（2,205 个细胞）和低丰度的 `1-InnateLymph`(64)/`2-Macrophage`(139)/`4-Myoid`(49) 在 bulk 里解不稳。
→ **主分析合并成 7 类**：`[SPG(8), Scytes(9), STids(10), Elongating(11), Sertoli(6), Leydig(5), Other(1,2,3,4,7)]`。
→ 敏感性分析保留 11 类。两版都跑，结论对比。

**NNLS 求解**
- 对每个 bulk 样本：`min ||S_marker @ p - y_marker||², p ≥ 0`，`scipy.optimize.nnls`。
- 每个样本的 p 归一化到和为 1。**不做无约束线性回归**（会出现负数比例）。
- 每个 bulk 数据集用**自己的基因交集**建 S 的子集（不同数据集 marker 数不同，正常）。
- 实测速度：2,000×9、80 样本 **0.02 s**。整件事不到 1 分钟。

**诊断（必须出，否则第二幕没法解释）**
1. 重建误差 `R² = 1 - ||y - S p||² / ||y - mean(y)||²`，逐样本画。
2. 残差与 S 的正交性检查（**必须检查的陷阱**：NNLS 的解按构造让残差与 S 的列正交，这会人为压低残差相似度——必须画出「残差与 S 各列的相关系数」并报告，若全部 ≈ 0 就是这个问题在起作用）。
3. 比例堆叠柱状图按样本排序，按数据集分组。
4. WT vs KO 的比例位移 `ΔP = mean(p_KO) - mean(p_WT)`，每个数据集一张。

**验收判据**
- `S_celltype.tsv` 是 **37,241 × 11**（不是 34,633 × 11——**行是基因**）。看一眼 `head`，第 1 列是基因名。
- `P_proportions.tsv` 每个样本的比例和 == 1（容差 1e-6），全部 ≥ 0。
- 重建 `R²` 中位数 **> 0.6**。低于 0.5 说明参考系与 bulk 不匹配，**这时第二幕的结论不可信**，必须报告而不是硬往下走。
- `ΔP` 方向符合生物学：KO 组 `Elongating`/`STids` 比例**下降**，`Sertoli`/`Leydig` 比例**上升**。**这是整个项目的正向对照**——如果这个方向都不对，去卷积就是错的。
- 峰值 RSS 打印出来，断言 **< 2 GB**。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| **内存爆 / R 报 `vector memory limit of 16.0 Gb reached`** | 说明脚本里有地方物化了矩阵。grep 脚本找 `read_csv`/`np.loadtxt`/`as.matrix`/`fread`，改成流式。**这是最可能翻车的点** |
| `MemoryError` / 机器卡死 | 降到「基因分块」方案：每块 1,000 个基因，`np.zeros((1000, 34633), dtype=np.float32)` = 139 MB/块，块间累加进 S。峰值 ~300 MB |
| 细胞数对不上（≠34,633） | 检查 `PerCellAttributes` 是否跳过了前 4 行注释。做断言 `len(labels) == 34633` |
| NNLS 解出负数或全 0 | 忘了加非负约束，或 markergene 太少。打印 `nnls` 的 `resid` |
| `R²` 很低（< 0.5） | 参考系物种/组织不匹配，或归一化尺度错了。试 `log2(CPM+1)` 空间下做 NNLS（而不是线性空间） |
| 某 bulk 数据集与参考系基因交集 < 300 | 该数据集平台太老（芯片）。**跳过它的去卷积**，只在第一幕用它 |
| Python 环境坏 | 用 `code/05_deconvolution.R`（`nnls` 包已装，实测 2,000×9 / 80 样本 0.03 s）。但 R 路径必须自己保证流式（`readLines` + `scan`），别用 `fread` |

---

### Step 6 · 第二幕：零模型检验（**核心**）

| 项 | 值 |
|---|---|
| 脚本 | `code/06_act2_null_model.R` |
| 输入 | `results/act1/`（相似度）、`results/deconv/P_proportions.tsv`（组成位移 ΔP）、`results/de/`（logFC） |
| 输出 | `results/act2/`（Mantel r 表、偏 Mantel 表、回归系数表、置换分布）、`act2/06.done` |
| 耗时 | Mantel 999 次置换 × 约 30 个数据集对 ≈ **2-5 min** |
| 峰值内存 | **< 500 MB**（全是小矩阵：8×8 距离、8×11 比例） |
| **可中断点** | ★ 每个检验（Mantel / 偏 Mantel / 回归）跑完后落盘一次 |

数学定义与陷阱处理见 **`docs/03-null-model.md`**——本文只写执行规格，不重复推导。执行上要落实这几条：

**要算的量**
1. `D_expr`：模型间 logFC 的差异矩阵（`1 - Spearman`）。
2. `D_comp`：模型间 ΔP 的差异矩阵。用 **Euclidean on CLR 变换后的 ΔP**（零模型规格 e 条：比例数据必须做变换，否则单纯形约束会让距离扭曲）。
3. **Mantel 检验**：`mantel(D_expr, D_comp)`，999+ 次置换。置换的是**样本标签**（模型身份的置换），不是矩阵元素——保持每个模型内部结构不变。
4. **偏 Mantel**：控制 `node`（三个阻滞节点），检验 ΔP 在扣掉节点后还有没有解释力。
5. **回归**：`dist_expr ~ dist_comp + node_same + platform_same + tissue_same`，用置换检验拿 p 值（不做正态假设）。
6. **方差分解**：`dist_expr` 的方差里，ΔP 贡献多少、node 贡献多少、platform 贡献多少。这是"组成能解释多少"的直接答案。

**验收判据**
- 置换分布直方图能画出来（999 个 null r 值），观测 r 在分布里的位置明确。
- 报告 `r` 的点估计 + 置换 p + 95% CI。
- **必须同时报告 node 和 platform 的贡献**——如果 platform 贡献比 ΔP 还大，那第一幕的"收敛"可能只是平台效应，这是个必须诚实报告的发现。
- 三组敏感性都跑：`S_global` vs `S_within_dataset`、7 类 vs 11 类、`unstable` 数据集剔除前后。结论跨版本一致才算稳。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| 模型数太少（< 6）做不了偏 Mantel | 先做简单 Mantel + 回归，偏 Mantel 等 Step 8 扩数据集后补 |
| 置换 p 全是 1.0 | 检查置换对象搞错了（别置换了 ΔP 的列标签，要置换模型对） |
| 距离矩阵不对称 | `mantel` 要求对称；`D = (D + t(D))/2` |
| ΔP 在 CLR 后出现 `-Inf` | 有比例恰好为 0。加伪计数（`+ 1e-6`）或用 `zCompositions` 的贝叶斯替换，**并在报告里说明** |

---

### Step 7 · 第三幕：残差化后重跑 DE 与跨模型比较

| 项 | 值 |
|---|---|
| 脚本 | `code/07_act3_residual.R` |
| 输入 | `data/interim/<GSE>_expr.rds`、`results/deconv/`（S、P） |
| 输出 | `results/act3/<GSE>_resid_de.tsv`（同 DE 表 schema）、`results/act3/resid_similarity.tsv`、`act3/07.done` |
| 耗时 | 8 数据集 × (组成预测 + 残差 + DE) ≈ **3-5 min** |
| 峰值内存 | **< 1.5 GB**（每数据集仍是 20,000 × 40 的小矩阵） |
| **可中断点** | ★ 每个数据集残差化并重跑 DE 后落盘一次 |

**做法（项目笔记 §5.3）**
1. 对每个样本：`组成预测表达 = S_marker @ p`（`p` 来自 Step 5）。
2. 观测（归一化到同一空间，`log2(CPM+1)`）减掉组成预测 → 残差矩阵（同尺寸，20,000 × N，6.4 MB，安全）。
3. 在**残差**上重跑 Step 3 的同一套 DE 流程 → `_resid_de.tsv`。
4. 在残差 logFC 上重跑 Step 4 的相似度 → 与第一幕对比。

**关键对照**：`results/act3/` 里必须同时有
- 第一幕的天真相似度（`results/act1/`）
- 残差化后的相似度
- **只扣 ΔP 的平均效应（不逐样本扣）** 的版本——用来区分"扣掉组成均值位移"与"逐样本校正"

**验收判据**
- 残差化后平均相似度**预期显著低于**第一幕（这是核心论点：相似度大部分是组成伪影）。
- 残差相似度**不应该是 0**——完全归零反而是可疑的（说明过度校正把信号也扣了，见 Step 8 校准）。
- **阳性对照必须做**：残差里 `Prm1`/`Tnp1` 这些精子细胞基因的表达应该**被扣掉**（因为它们是组成驱动的）；如果残差里这些还是最显著的下调基因，说明残差化没生效。
- 三幕的三张相似度矩阵并排存好，Step 9 画成大图。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| 残差相似度 ≈ 0 或为负 | 过度校正。看 Step 8 校准结果；考虑用「平均 ΔP 校正」而非「逐样本校正」 |
| 残差里残差方差异常小 | 检查 `S @ p` 的尺度有没有比观测大（归一化空间不一致，最常见是 log 与线性搞混） |
| 某数据集没有去卷积结果（Step 5 跳过了） | 该数据集在第三幕里跳过，并在报告里列出 |
| 残差出现负值 | 正常（log 空间相减），但 DESeq2 吃不了负数 → 残差必须用 `limma` 跑，**不能用 DESeq2** |

---

### Step 8 · 校准与敏感性分析（伪 bulk 对照）

| 项 | 值 |
|---|---|
| 脚本 | `code/08_calibration.py` |
| 输入 | `results/deconv/S_celltype.tsv`、`results/deconv/P_proportions.tsv`、GSE112393 原始 DGE（复用 Step 5 的流式读取函数） |
| 输出 | `results/calib/`（回收误差、注入-恢复曲线、敏感性网格）、`calib/08.done` |
| 耗时 | 合成 100 个伪 bulk（每个 = 对 34,633 细胞的抽样 + 加权求和）≈ **5-10 min** |
| 峰值内存 | **< 500 MB**（伪 bulk 是 20,000 × 1 的向量，逐个生成） |
| **可中断点** | ★ 每完成一组（例如「每组 10 个伪 bulk」）落盘一次 |

**三件事（对应零模型规格里的阳性对照设计）**

1. **比例回收**：从参考系里按已知比例 `p_true` 合成伪 bulk → 跑 NNLS → 比较 `p_hat` 与 `p_true`。成功判据：中位绝对误差 < 0.05，且每个类型的回收误差可接受。
   - **注意内存**：合成伪 bulk 时不要 `X @ p`（X 是全矩阵）。用流式：逐基因读一行，`np.dot(row, p)` 直接得该基因的伪 bulk 值。峰值仍然是单行 277 KB。
2. **注入已知效应**：在选定的 marker 基因上注入已知 logFC（例如 1.0、2.0 倍），看残差化后能不能回收。这直接回答"我们会不会把真信号也扣掉了"。
3. **敏感性网格**：
   - 参考系：`S_global` vs `S_within_dataset`
   - 粒度：7 类 vs 11 类
   - 归一化：CPM vs TPM vs 中位数
   - 空间：线性 vs log
   - 数据集子集：全部 vs 仅 n≥6 vs 仅 RNA-seq
   每格跑一遍完整的第一/二/三幕，记录结论是否翻转。

**验收判据**
- 比例回收：`p_hat` vs `p_true` 的散点要贴对角线，报告 R² 和最大绝对误差。
- **注入-恢复斜率接近 1**。如果斜率明显 < 1，说明残差化吃掉了真信号，第三幕的所有"剩余信号"结论都要打折。
- 敏感性网格里，主结论（"组成解释大部分相似度"）在**≥ 80% 的格子里方向一致**。不一致的组合必须在报告里点名。
- 校准**失败是被允许的结果**（零模型规格的结果 C），但必须显式写出来，不能跳过。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| 比例回收 R² < 0.7 | 参考系粒度与 bulk 不匹配。合并细胞类型后重试 |
| 注入的信号回收不到 | 残差化过强。改用零模型规格建议的「平均 ΔP 校正」 |
| 伪 bulk 合成时内存爆 | 一定是不小心 `X @ p` 了。改成逐行 `dot` |
| 敏感性网格跑不完（时间不够） | 先跑 3 个最关键的轴（参考系、粒度、归一化），其余记 TODO |

---

### Step 9 · 出图与出表

| 项 | 值 |
|---|---|
| 脚本 | `code/09_figures_tables.R` |
| 输入 | `results/` 下全部中间结果 |
| 输出 | `figures/fig01..figNN.*`、`results/tables/table01..NN.tsv`、`09.done` |
| 耗时 | **5-10 min**（纯画图） |
| 峰值内存 | **< 1 GB** |
| **可中断点** | ★ 每张图独立函数 + 独立落盘，跑完一张存一张 |

**图清单**

| 编号 | 画什么 | 用什么数据 | 图类型 |
|---|---|---|---|
| fig01 | 数据集总览：三个节点 × 样本量 × 平台 | `data/interim/*_meta.tsv` | 分面条形 / 点图 |
| fig02 | 第一幕：跨模型相似度热图 + 聚类树 | `results/act1/` | 热图（`pheatmap`） |
| fig03 | 第一幕：top-N 基因重叠（Jaccard vs N） | `results/act1/` | 折线 |
| fig04 | 去卷积诊断：重建 R² 分布 | `results/deconv/` | 直方图 / 箱线（按数据集分面） |
| fig05 | 去卷积诊断：残差 vs S 各列的相关（正交性陷阱） | `results/deconv/` | 热图 |
| fig06 | 细胞比例堆叠柱状图（按样本，按数据集分面） | `results/deconv/P_proportions.tsv` | 堆叠柱状 |
| fig07 | ΔP 森林图：每个数据集的 WT→KO 比例位移 | `results/deconv/` | 森林图 / 点线图 |
| fig08 | 第二幕：Mantel 置换分布 + 观测值 | `results/act2/` | 直方图 + 竖线 |
| fig09 | 第二幕：方差分解（ΔP / node / platform） | `results/act2/` | 堆叠条形 / 韦恩式 |
| fig10 | **三幕并排**：天真 / 组成校正 / 残差的相似度矩阵三联图 | `act1/` + `act3/` | 三个热图拼版（`patchwork`） |
| fig11 | 残差 DE 的火山图（每个数据集一张，分面） | `results/act3/` | 火山图 |
| fig12 | 校准：`p_true` vs `p_hat` 散点 | `results/calib/` | 散点 + 对角参考线 |
| fig13 | 校准：注入-恢复斜率 | `results/calib/` | 散点 + 拟合线 |
| fig14 | 敏感性网格：结论是否翻转 | `results/calib/` | 热图或点阵 |
| fig15 | 组成驱动基因 vs 残差剩余基因的对比（Prm1/Sox9 等标记基因） | `de/` + `act3/` | 箱线 / 小提琴 |

**表清单**

| 编号 | 内容 | 来源 |
|---|---|---|
| table01 | 数据集清单（GSE、模型、节点、n、平台、supp 类型） | `data/dataset_selection.csv` + `_meta.tsv` |
| table02 | 每数据集 DE 汇总（方法、显著基因数、线粒体占比前后） | `results/de/_index.tsv` |
| table03 | 共同基因数与各相似度指标 | `results/act1/` |
| table04 | 细胞类型特征谱（S）的 marker 基因 | `results/deconv/markers.tsv` |
| table05 | 各样本细胞比例（P） | `results/deconv/P_proportions.tsv` |
| table06 | Mantel / 偏 Mantel / 回归结果 | `results/act2/` |
| table07 | 方差分解表 | `results/act2/` |
| table08 | 残差化前后相似度对比 | `act1/` + `act3/` |
| table09 | 校准结果（回收误差、注入-恢复） | `results/calib/` |
| table10 | 敏感性分析结论对照 | `results/calib/` |

**验收判据**
- 15 张图全部生成，每张都是矢量或 ≥ 300 dpi（投稿用）。
- 每张图的**数据来源都能追到某个 `.tsv`**（在图脚本里写注释标注来源路径）。
- table01-table10 全部是 TSV，列名自解释。
- **图里不要写"显著""重要"这类词**——措辞纪律见项目笔记 §3.1「我们检验这类分析的前提假设」。

**如果这步失败怎么办**
| 症状 | 处理 |
|---|---|
| `ComplexHeatmap` 缺失 | 已实测缺失，**用 `pheatmap`**（已装）。不要尝试装 |
| 中文乱码 | 图里用英文标签（投稿本来也要英文） |
| 某张图数据缺 | 该图函数 try-catch 跳过并记日志，不要因为一张图失败中断全部 |

---

## 5. 总耗时汇总

### 5.1 机器计算时间（不含下载，不含人思考）

| Step | 内容 | 耗时 | 峰值内存 | 可中断点 |
|---|---|---|---|---|
| 0 | 环境自检 | < 0.5 min | < 0.5 GB | — |
| 1 | 下载（仅参考系） | **11 min** | ~0.1 GB | ★ 每文件 |
| 1 | 下载（+8 个 bulk，视清单） | +5-30 min | ~0.1 GB | ★ 每文件 |
| 2 | 读入 + 基因统一 | 2-4 min | < 1.5 GB | ★ 每数据集 |
| 3 | 每数据集 DE | 1-2 min | < 1.2 GB | ★ 每数据集 |
| 4 | 第一幕相似度 | < 1 min | < 0.5 GB | 跑完即止 |
| **5** | **去卷积（流式）** | **~5 min** | **< 0.5 GB**（实测 47.7 MB） | ★★ 三处 |
| 6 | 第二幕零模型 | 2-5 min | < 0.5 GB | ★ 每检验 |
| 7 | 第三幕残差 DE | 3-5 min | < 1.5 GB | ★ 每数据集 |
| 8 | 校准 + 敏感性 | 5-10 min | < 0.5 GB | ★ 每 10 个伪 bulk |
| 9 | 出图出表 | 5-10 min | < 1 GB | ★ 每张图 |
| | **合计（计算）** | **约 40-60 min** | **全局峰值 < 1.5 GB** | |

> 除 Step 5 的一次性 5 分钟外，**没有任何一步超过 10 分钟**。这意味着每一步都能塞进「每周 2 小时」的一个碎片里——**一个会话可以跑完 3-5 个 step**（含人核对元数据的时间）。

### 5.2 人类时间（真正的时间瓶颈）

| 活动 | 估计 |
|---|---|
| 手工核对每个数据集的 sample_titles → WT/KO 分组 | **每个数据集 5-15 min，8 个 = 1-2 h**（这是最耗时也最容易出错的环节） |
| 确认 supp 文件到底是计数矩阵还是作者 DEG | 每个 5 min，8 个 = 40 min |
| 读 `docs/03-null-model.md`、理解 CLR/偏 Mantel | 1-2 h |
| 看图、判断结论是否成立、写报告 | 4-8 h |

**机器只占 1 小时，人占 8-15 小时。** 排期要按人算，不是按机器算。

### 5.3 建议的会话拆分（每周 2 小时）

| 会话 | 内容 | 机器时间 | 人的时间 | 断点 |
|---|---|---|---|---|
| S1 | Step 0 + 启动 Step 1 下载（后台） | 11+ min | 20 min | 下载在后台跑，下次收 |
| S2 | Step 1 收尾 + 核对 supp 类型 + Step 2 | 5 min | 60 min | 元数据核对 |
| S3 | Step 2 收尾 + Step 3 + Step 4 | 5 min | 40 min | — |
| S4 | **Step 5 去卷积** | 5 min | 40 min | 诊断图判读 |
| S5 | Step 6 + Step 7 | 5 min | 60 min | — |
| S6 | Step 8 校准 | 10 min | 30 min | — |
| S7 | Step 9 出图出表 | 10 min | 60 min | — |

**7 个会话 × 2 小时 ≈ 3.5 周**。加上写报告和读文献，与项目笔记 §6 的「4-6 周」一致。

---

## 6. 内存风险排序与降级方案

### 🔴 风险 1（最高）：Step 5 去卷积——参考系稠密读入

| 项 | 内容 |
|---|---|
| **为什么危险** | 37,241 × 34,633 × 8 = **10.32 GB**。R 自己的向量堆上限是 **16 GiB**（实测 `mem.maxVSize()` = 16384 Mb）。实测：分配成功，但**任何一次复制都直接报错** `vector memory limit of 16.0 Gb reached`。所以不是"慢"，是"必死" |
| **触发条件** | 脚本里出现 `read_csv`/`np.loadtxt`/`fread`/`as.matrix` 任何一个读全矩阵的调用 |
| **主方案** | 流式聚合，**不物化矩阵**。实测峰值 **47.7 MB** |
| **降级 A** | 「基因分块」：每块 1,000 基因 × 34,633 细胞 × float32 = **139 MB/块**，块间累加进 S。峰值 ~300 MB。**注意用 float32，float64 会翻倍到 277 MB/块** |
| **降级 B** | 只读 2,000 个 marker 基因：2,000 × 34,633 × float32 = **277 MB**。但需要事先知道 marker 是哪 2,000 个——所以只在"已知 marker 列表"时才用 |
| **降级 C** | 用 `GSE112393_12GermCellClusters_AllGeneExp.xlsx`（11.4 MB，作者已聚合到 12 个生殖细胞簇）替代原始 DGE。代价：只有生殖细胞、没有体细胞，无法解出 Sertoli/Leydig——**这会让第二幕的核心论证残缺**，只在救急时用 |
| **验证方式** | 脚本里加 `resource.getrusage(RUSAGE_SELF).ru_maxrss` 打印，断言 < 2 GB |

### 🟠 风险 2：下载遇到 496 GB 的陷阱文件

| 项 | 内容 |
|---|---|
| **为什么危险** | `GSE244681_RAW.tar` 实测 **496,539.8 MB = 496 GB**。磁盘只有 82 GB。按 126 KB/s 要 **64,140 分钟（44 天）** |
| **主方案** | 下载清单**先 probe 再人工过目**（Step 1 第 1-2 步）。任何 `RAW.tar` > 300 MB 默认不下 |
| **降级** | 该数据集放弃，用同节点的替代数据集补位 |
| **验证** | Step 1 里加断言：清单里任何 `size_bytes > 3×10^8` 且 `needed_for` 不含 `required` 的行直接拒绝 |

### 🟢 风险 3（已消除）：`fetch_geo.py` 下载可靠性

> **状态：已修复，不再是风险。** 早期版本确实有两个缺陷——分块攒在内存 `parts` 字典里、中断即全丢；以及缺少"目标文件已完整"的守卫，导致陈旧的 `.parts` 目录会把一份好文件重建成截断版本（参考矩阵就是这样被弄坏的）。**这两个缺陷在当前版本（211 行）中都已修复，并已实测验证。**

| 项 | 内容 |
|---|---|
| **当前防护 1** | `scripts/fetch_geo.py:139-149` 每块到达后立刻写自己的 `.tmp` 再 `os.replace` 成正式块 → 中断只丢当前块，重跑续传 |
| **当前防护 2** | `scripts/fetch_geo.py:121-137` 每轮扫描已有块并跳过尺寸正确的 → 真正的块级续传 |
| **当前防护 3** | `scripts/fetch_geo.py:98-103` 目标文件已存在、尺寸正确、且 `gzip_ok()` 通过 → 打印 `already complete ... -- skipping` 直接返回，**绝不动它**。这正是当初损坏参考矩阵的那个 bug 的修复 |
| **当前防护 4** | `scripts/fetch_geo.py:184-187` 组装后对 `.gz` 做 `gzip_ok()` 全量 CRC 校验，不通过则**拒绝安装**（`refusing to install it`），不会留下坏成品 |
| **当前防护 5** | `scripts/fetch_geo.py:114-117` 启动时清理残留的 `.tmp` 草稿；`:194-196` 成功后 `rmtree` 掉 `.parts` |
| **当前防护 6** | `scripts/fetch_geo.py:60-62` 校验 `len(data) == want`，503 的短正文被判 `short read` 并重试，不会污染分块 |
| **内存** | 上界 = `streams × chunk` = 6 × 3 MB = **18 MB**，与文件总大小无关 |
| **残余注意** | `gzip_ok()` 会对整个文件做一次完整解压（3.7 GB payload 约 1-2 min），这是刻意的成本——它正是能抓住截断的那个检查 |

**实测验证**（4 块 × 400 KB，`--chunk 400000`）：完整下载与原始 `cmp` 字节一致；删块重跑只补 1 块且成品一致；把块截成 0 字节后重跑被识别为缺失并重下；HEAD 404 时不创建垃圾文件。
→ **Step 1 不需要为它写任何降级方案。** 唯一要做的运维动作是下载完成后确认 `.parts` 已被清掉（正常情况下脚本自己会删）。

### 🟢 风险 4：DESeq2 在样本多的数据集上

| 项 | 内容 |
|---|---|
| **实测** | 20,000 × 12 峰值 **0.93 GB**。外推到 n=40 约 2-3 GB |
| **判断** | 安全，但 16 GB 下同时开多个 R 进程会叠加 |
| **降级** | 串行跑；或改用 limma-voom（0.6 s，内存更低） |
| **注意** | **绝不要**并行跑 8 个数据集的 DESeq2（8 × 2 GB = 16 GB 直接满） |

### 全局原则

- 任何时刻**所有大对象之和 ≤ 6 GB**，留 10 GB 给 R/OS 的临时副本。
- 用 `.venv/bin/python`，**不要**用 `.venvtest`（只有 numpy）或 conda base（没有 numpy）。
- 大对象用完立刻 `rm(); gc()`。

---

## 7. 最容易被卡住的环节（按卡住概率排序）

1. **🔴 Step 5 的稠密读入（概率最高）**。三个月后手滑写一句 `fread(dge)` 就完蛋，而错误信息（`vector memory limit of 16.0 Gb reached`）看着像"内存不够"，容易误判成"加内存"或"换机器"，而不是"改读法"。**在 `code/05_deconvolution.py` 文件头写一句大字注释：本文件禁止出现 loadtxt/read_csv/fread。**

2. **🔴 元数据分组核对（最耗人时）**。`data/inventory.csv` 的 `n_ko`/`n_wt` 已知有错（GSE235777 记 3/6）。如果 WT/KO 搞反，整条流水线跑得很顺，结论完全反向，而且**不会报错**。唯一防线是 Step 3 的验收判据：「Prm1 应下调、Sox9 应上调」。**这个断言必须写成硬性检查，不能靠人看。**

3. **🟠 NCBI 503 限流（续传已修复，限流仍在）**。503 会成片持续（实测连续 6 次、退避到 90 s 仍失败），且返回的是 HTML 错误页。`fetch_geo.py` 已能识别并重试、且可断点续传（§2.6 风险 3），所以**下载不会白费**；但**它可能连续几小时下不完**。会话制下（每周 2 小时）**建议 S1 就把下载挂后台（`nohup`），别占着会话**，下次开工先看 `.parts` 是否已清、成品是否 `gzip -t` 通过。

4. **🟠 `RAW.tar` 里的文件类型与预期不符**。GEO 的 supp 里经常混着 ChIP-seq/Hi-C 文件（GSE244681、GSE242235、GSE252121 都这样）。`data/inventory.csv` 的 `matrix` 字段标了 `有（补充文件含矩阵）`，但**不等于**那个矩阵就是我们要的 bulk 计数矩阵。**下完先 `tar -tf` 看一眼成员名**。

5. **🟡 GSE259417 只有作者 DEG**。实测 supp 只有 `GSE259417_P16_DEG.csv.gz`（2.4 MB），是作者已经算好的差异表。它能进第一幕（有 logFC），**但不能进去卷积和残差化**（没有原始计数就无法建模组成）。**必须在清单里显式标记，否则 Step 5 会拿不到它的表达矩阵而报错。**

6. **🟡 矩阵朝向/对数空间搞混**。去卷积要求线性强度空间，DE 要求 log 空间。同一份数据在不同 step 用不同变换，是隐蔽的 bug 源。**在第 2 步的元数据里显式记录 `quant`（counts/FPKM/TPM/CPM）和 `log_state`（raw/log2），后面每步都读这个字段决定变换。**

7. **🟡 R 版本升级导致 `mem.maxVSize` 变化**。Step 0 里断言它 == 16384；如果变了，§3 的 6 GB 预算要重算。

---

## 8. 一页速查（贴在显示器上）

```
每次开工：
  bash code/00_check_env.sh                 # 30 s，确认环境没坏
  ls results/*/*.done                       # 看跑到哪了

内存铁律：
  禁：读入 37241 × 34633 稠密矩阵（10.32 GB，R 上限 16 GiB，复制即死）
  行：流式聚合到 S（37241 × 11 = 3.3 MB），实测峰值 47.7 MB

工具：
  .venv/bin/python        ← 分析用（有 numpy/pandas/scipy）
  Rscript + .libPaths(c(".Rlib", .libPaths()))
  python3 scripts/fetch_geo.py <url> <dest>   ← 6 并发，3 MB/块，可续传；重跑即续，别删 .parts

进度哨兵：
  data/interim/02.done   results/de/03.done    results/act1/04.done
  results/deconv/05.done results/act2/06.done  results/act3/07.done
  results/calib/08.done  figures/09.done

最可能翻车的两句：
  Step 5 里 "先把矩阵读进来" → 必死，改流式
  Step 3 后 "Prm1 应该下调"  → 没下调就是 WT/KO 反了，别往下跑
```
