# 生精阻滞转录组共性分析 —— 本机可执行方案

> 适用机器：Apple M4 / 10 核 / **16 GB RAM** / 可用磁盘 ~80 GB / macOS aarch64
> 本文件替代 `notes/原计划书.md`。所有"能跑"的结论都经过本机实测，实测证据见 `docs/01-environment.md` … `docs/04-pipeline-runbook.md`。

---

## 一、先给结论

**这份文件夹里有两个互相矛盾的文档，方案必须建立在正确的那一个上。**

| 文档 | 状态 | 判断 |
|---|---|---|
| `notes/原计划书.md`（214 行） | ❌ **数据层面不成立** | 「≥3 个 miRNA KO 共有差异 = 核心基因」需要 ≥3 个 miRNA KO 转录组。实测**只有 miR-202 一个**。n=1 做不了交集，也做不了相似度矩阵。 |
| `notes/思路整理.md`（219 行） | ✅ **方向正确，且地基已验证** | 本项目的早期笔记已推翻原计划书，把问题从「找共性」改成「**这些共性是真的吗？多大程度上只是细胞组成变化的副产物？**」 |

**本次方案采纳 `notes/思路整理.md` 的科学问题**，并已经用真实数据把它的地基跑通了：

> **正向对照已通过**：用 GSE236388 的 XAP5 cKO 睾丸臂做端到端去卷积，KO 组 germ-cell 比例下降 **0.194**、Somatic 比例上升 **0.194**，重建 R² = 0.680–0.742（阈值 0.6）。
> 说明参考系能回收一个**生物学上方向正确**的组成位移。第二幕（核心）有可用地基，不是纸上谈兵。
> 复现命令：`code/06_validate_deconvolution.py`，原始输出见 `docs/02-datasets.md`。

**一句话方案**：用 8 个**全睾丸**生精阻滞 bulk 数据集 + GSE112393 单细胞参考图谱，先按常规跑一遍跨模型收敛（第一幕），再构造**只含细胞组成**的零模型去量化"组成能解释多少相似度"（第二幕），最后在残差里找真正剩下的信号（第三幕）。
**机器只占约 1 小时；主要成本是 8–15 小时的人时**（元数据核对 + 看图写报告）。

---

## 二、为什么原计划书不成立（实测证据）

| 原计划书的假设 | 实测结果 |
|---|---|
| 有 6 类 miRNA KO 小鼠睾丸转录组可用 | 逐个核查后**只有 miR-202**（GSE174731 / GSE126936 / GSE99160） |
| miR-34b/c+miR-449、miR-506、XmiRs 可用 | **查不到**对应数据集 |
| miR-17-92（GSE225675/GSE225676）可用 | 是 **XY 完全性别反转**模型，与生精阻滞不可比 |
| 可做 Jaccard 相似度矩阵 + 聚类 | n=1 **无法**做交集/相似度矩阵 |
| 参考文献第 3、5 条 | **重复**（同一篇 Wu et al. 2014 PLoS Genet 列了两次） |

**技术路线上还有两个隐藏问题**（原计划书未察觉）：

1. **RNA-seq + 芯片混用会引入平台伪影**。「跨模型共性」很可能一部分来自平台批次效应，而非生物学。
2. **细胞组成混杂**。睾丸 bulk 里所谓"共识信号"（`Prm1/Prm2/Tnp1/Tnp2` 下调、`Sox9/Wt1/Gma4/Amh/Cyp11a1` 上调）在生物学上就是「**失去了晚期生精细胞**」这一件事的另一种说法 —— 精原细胞阻滞模型里这些基因本来就会变。这是本方案要检验的核心，不是要回避的偏差。

---

## 三、三幕结构（方案骨架）

```
第一幕  天真版：跨模型 logFC 收敛
        └─ 复现"不同阻滞模型看起来很像"这个文献常见观察
第二幕  核心：只含细胞组成的零模型          ← 本方案的价值所在
        └─ 问：只用细胞组成能解释多少跨模型相似度？
第三幕  残差里还有什么
        └─ 扣掉组成预测后，真正剩下的信号
```

### 第二幕的数学定义（这是整个项目最容易被打穿的一环）

**零模型预测的 logFC**（组成位移经参考谱推进到表达空间、且必须带 log）：

$$\hat c_m(g)=\log\!\big((S\,\bar p^{KO}_m)(g)+1\big)-\log\!\big((S\,\bar p^{WT}_m)(g)+1\big)$$

含义：**如果每个基因的表达量本身毫无变化、只有细胞组成变了**，模型 m 会「看起来」有哪些基因上调/下调。

**三个用模拟实测出来的坑（写进方案以避免踩）**：

1. 🔴 **最自然的估计量是错的。** 直接用「logFC 投影到 S 的列空间」作组成预测，**即使给完美参考系也只能报 0.585（真值 1.0），系统性低估约 40%** —— 因为 log 是凹函数、组成效应不落在 S 的列空间里。改用上式（经 S 推进的组平均比例预测）后回到 **0.992**。
2. 🔴 **朴素 Mantel 检验无信息。** 模拟实测：在「纯伪影」和「有真信号」这两种**相反**情形下，Mantel 给出**同一个 p = 0.158**。已弃用，改用**以模型对为单位的 MRM 回归**。
   - 三条理由：① Legendre 2015 已系统批评 Mantel 并推荐 MRM；② Mantel 操作样本间距离矩阵，而 ΔP 是每个模型**内部**的 KO-vs-WT 对比量；③ 模拟实测对核心问题无判别力。
3. 🔴 **置换的单位必须是「模型」而不是基因或样本**，且必须用**模型错配**（把模型 i 的 logFC 与模型 j 的 ĉ 配对）。实测基因打乱法的零分布 sd = 0.007 **过窄**，会大量假阳性；模型错配法 sd = 0.129 才是正确的。
   - 置换次数：M=6 穷举 720、M=7 穷举 5040、M=8 优先穷举 40320。

**主检验**：

$$y_{ij}=\mathrm{cov}(\mathrm{lfc}_i,\mathrm{lfc}_j)\sim \mathrm{cov}(\hat c_i,\hat c_j)+\mathbb 1[\text{node}_i=\text{node}_j]+\cdots$$

### 🔴 必须诚实报告的短板：检验力

| M（模型数） | a=0（纯伪影） | a=0.20 | a=0.40 |
|---|---|---|---|
| 6 | 崩溃（1.59, **sd 19.9**） | 0.293 ± 0.158 | 0.103 ± 0.061 |
| 12 | — | sd 降到 0.03 | — |

⇒ **compShare 一律报「点估计 + 置换 95% 区间」，绝不报两位小数。**
**M=6–8 能回答「组成是不是主要解释」，不能回答「37% 还是 44%」。** 这句话必须写进论文 Limitations，否则会被审稿人抓住。

### 措辞纪律（`notes/思路整理.md` 已强调，方案沿用）

❌ 不写「我们指出前人错了」
✅ 写「我们**检验这类分析的前提假设**」

最近先例 **Zhang J et al., "Microarray meta‐analysis identifies candidate genes for human spermatogenic arrest", *Andrologia* 2019;51(8):e13301**（PMID 31037746，DOI 10.1111/and.13301；Crossref 标题已实测核对）—— 标题本身就说明它做的正是「天真版」：178 例人睾丸 / 7 套芯片 / 按 Johnsen 评分分 4 型做 meta 分析取交集。
**注意**：该文全文在 Wiley 付费墙后，**未取得**（PMC 无副本、Europe PMC `isOpenAccess:N`）。摘要中**没有**任何关于细胞组成/去卷积的讨论，但**正文不可得故无法确认** —— **不要声称该文"忽略了组成"**。另：PubMed 记录首作者为 **Zhang J**，与 `notes/思路整理.md` 写的 "Kui F" 不符。

---

## 四、数据：CORE-8（全部实测为全睾丸）

**硬筛子**：去卷积**只能用于全睾丸样本**。分选细胞群体的组成由分选决定、不是测量出来的，其 P 不是生物学估计，残差化无意义。`code/05_deconvolution.R` 里有硬性 `stop()` 守卫，非 `whole_testis` 直接拒绝。

| # | GSE | 模型 | 节点 | n (KO/WT) | 年龄 | 平台 | 去卷积 |
|---|---|---|---|---|---|---|---|
| 1 | **GSE38096** | Sin3a 生殖细胞特异 cKO (Vasa-cre) | 精原细胞/起始 | 6 (2/2+2HET) | P0 | GPL13912 Agilent 4x44K | ⚠️ 需归一化 |
| 2 | **GSE42601** | Dazap1-Fn（MI 前阻滞） | 精原细胞/起始 | 6 (3/3) | P21 | GPL1261 Mouse430_2 | ⚠️ 需 RMA |
| 3 | **GSE294209** | Rnu11 cKO | 减数分裂前期 | 6 (3/3) | P14 | GPL24247 NovaSeq | ⚠️ TPM 非 counts |
| 4 | **GSE236388** | XAP5 cKO **仅 D16 睾丸臂** | 减数分裂前期 | 6 (3/3) | D16 | GPL24247 NovaSeq | ✅ **首选** |
| 5 | **GSE101511** | Hnrnpgt cKO | 减数分裂前期 | 6 (3/3) | 18 dpp | GPL19057 NextSeq | ⚠️ 转录本→基因 |
| 6 | **GSE283904** | N6amt1-SKO | 减数分裂前期 | 8 (4/4) | 成年 | GPL24247 NovaSeq | ✅ **首选** |
| 7 | **GSE68283** | Rfx2-/- | 减数分裂后/精子形成 | 9 (4/5) | P29–30 | GPL17021 HiSeq | ✅ **首选** |
| 8 | **GSE68284** | Rfx2-/-（同论文 P21） | 减数分裂后/精子形成 | 11 (6/5) | P21 | GPL17021 HiSeq | ✅ **首选** |

节点覆盖 **2 / 4 / 2**；**8/8 全部 `tissue_class=whole_testis`**；可去卷积 5 个 yes + 3 个 yes-with-caveat、**0 个 no**。逐条 GEO 原文证据在 `data/dataset_selection.csv` 的 `whole_testis_evidence` 列。

### ⚠️ 库存有 15 处标错（`data/inventory.csv` 不可直接信）

标"全睾丸"但**不是**全睾丸的：
- **GSE323391** = 睾丸**类器官**（体外）
- **GSE237030** = 10X scRNA + ES 诱导体细胞
- **GSE125495** = **雌性**性腺体细胞
- **GSE310340** = 系列主题确是 JQ1 抑制减数分裂前期，**但该号下的 20 个 GSM 实测为 Kidney×10 + Liver×10**（`Sample_source_name_ch1` 计数；应为毒性臂子系列）→ 此号本身不含睾丸样本，不可用
- **GSE23119** = **分选**精原细胞（BSA 梯度，纯度 >90%）
- **GSE236388** 的 n=12 是**混装**：6 个 D16 睾丸 + 6 个 D56 **附睾尾精子** → **只有 6 个可用**
- **GSE221433** 全睾丸正确，但模型是**西方饮食**不是遗传阻滞（作对照，不作 case）

标"待查"实为分选：GSE158471 / GSE253154 / GSE288152 / GSE324232
`supp_state=matrix` 标错：GSE259417（只有作者 DEG）、GSE114042（作者 DEG 表）、GSE57745（转座子家族计数，非基因水平）
计数错误：GSE235777 实为 3KO/3WT（库存写 3/6）；GSE101511 实为 3/3（库存为空）

**GSE259417**：只有作者 DESeq2 DEG 表（`KO_rep2` 还缺失）→ **能算 logFC ✅、不能去卷积 ❌**。

### 补充：EXTENDED-4（备用）
GSE213442 (Kpna2)、GSE214201 (Iqch，矩阵 4WT/3KO)、GSE235777 (Msl3)、GSE221433 (Tas1r3，饮食对照)

### 下载量（不是瓶颈）
| 口径 | 字节 | @126 KB/s |
|---|---|---|
| CORE-8 supp | 121,550,497 | 0.26 h |
| GSE112393 参考系 | 82,756,793 | 0.18 h |
| **合计** | **204,307,290** | **≈ 27 分钟** |

⚠️ 反例：`GSE244681_RAW.tar` 单个 **462 GB**（按 126 KB/s 要 **44 天**）—— 这类必须放弃。

---

## 五、为什么这台 16 GB 机器能跑（实测）

### 5.1 🔴 头号内存风险：参考系稠密读入必死

GSE112393 真实尺寸（**从 gzip 文件头读出，不是估的**）：**37,241 genes × 34,633 cells**。

| 方案 | 内存 | 结果 |
|---|---|---|
| 稠密读入（float64） | **10.32 GB** | ❌ `t()`/`scale()`/子集都会再造一份 → 20.6 GB |
| 稠密读入（float32） | 5.16 GB | ⚠️ 勉强，但没有余量 |
| 基因行分块 1000 行 | 139 MB/块 | ⚠️ 可行但没必要 |
| **流式聚合到 S（采用）** | **实测峰值 47.7 MB** | ✅ |

**实测的失败模式**（精确报错，供识别）：
```r
G <- 37241; N <- 34633
m <- matrix(0, nrow=G, ncol=N)   # ✅ 成功，object.size 报 9.6 Gb
m2 <- m * 2                      # ❌ cannot allocate vector of size 9.6 Gb
                                 #    vector memory limit of 16.0 Gb reached, see mem.maxVSize()
```
⚠️ **这条报错看着像"内存不够"，很容易被误判成"该换机器"** —— 其实只是方法错。换机器不解决，换方法才解决。

- **R 自己的硬墙**：`mem.maxVSize()` = **16,384 Mb = 16 GiB**（R 的向量堆上限，不是操作系统的）
- **内存安全边界**：参考系基因子集 5,000 → 1.39 GB ⚠️ / 子集 2,000 → 554 MB ✅ / **流式产出的 S（37,241×11）= 3.3 MB ✅ 主路径** / 单 bulk 数据集（20,000×40）= 6.4 MB / DESeq2 全开销实测峰值 **0.93 GB**
- ⚠️ **绝不要并行跑 8 个数据集的 DESeq2**（8 × 2 GB = 16 GB 直接满）

### 5.2 环境（conda 已弃用，走 venv + 项目内 .Rlib）

- **系统 R library 不可写**（`Operation not permitted`）→ 必须用项目内 `.Rlib`：`.libPaths(c(".Rlib", .libPaths()))`。**已装 171 个包 / 934 MB**，21 个关键包 `library()` 全部 PASS。
- **关键纠正**：不是"完全没有预编译二进制"。实测 CRAN 与 Bioc 3.23 **都命中 arm64 预编译二进制**（`bin/macosx/sonoma-arm64/contrib/4.6/`），21 个包**零编译痕迹**，**不依赖 clang**，安装远快于预期。
- `MuSiC` / `bisque` 装不上（GitHub 6 次重试全败）→ **用 `nnls` + `quadprog` 自实现 NNLS**，配合已装的 `EpiDISH` 作备选。
- **Python 3.14.7 上 scanpy 有 wheel**；项目内 `.venv` 已有 numpy/pandas/scipy/sklearn/openpyxl/matplotlib。
- ⚠️ **conda 并非完全不可用**（env-verifier 做了 A/B 对照推翻了我的初判）：重定向 `XDG_CACHE_HOME` + `CONDA_PKGS_DIRS` + `CONDA_REGISTER_ENVS=false` 后 conda 正常工作。但失败的根因只是写 `~/.conda/environments.txt` 这一个注册文件 → **本方案仍不用 conda**，理由从「不可用」改为**工程取舍**（venv + .Rlib 已可用、重定向脆弱）。workaround 见 `docs/01-environment.md`。
- ⚠️ **PyPI 直连被 connection reset，必须用清华镜像**（约 2.6 MB/s）—— 已内置进 `envs/setup.sh`。
- ⚠️ **matplotlib 需要 `MPLCONFIGDIR="$PWD/.mplcache"`**（`~/.matplotlib` 不可写）。

### 5.3 下载（GEO 限速，必须并发）

NCBI GEO FTP **单流仅 15–35 KB/s**。实测并发聚合：1 流 21 KB/s / 4 流 61 KB/s / **8 流 126 KB/s（最优）** / 16 流 81 KB/s。
`aria2c`/`wget`/`axel`/`lftp` **都没装** → 用 `scripts/fetch_geo.py`（stdlib，8 并发 range，3 MB/块，**块级可续传**）。

⚠️ **两个已经踩过并修好的坑**：
1. macOS/BSD 的 `seq` 生成大范围分块会把它解释成科学计数法 → `value too great for base (error token is "2e")`，静默产生 0 字节文件。**改用 Python 生成分块。**
2. 早期版本把所有块攒在内存里最后才写盘 → 一次瞬时 `HTTP 503` 让 20 分钟下载全丢。**已改为每块立刻落盘 + 多轮重试 + 指数退避**。并且：拼装成功后必须 `rmtree(.parts)`，否则残留的残缺 part 会在下次运行时**覆盖完好的文件**（真实发生过，把 81,261,540 B 的参考系写成了 79,322,404 B）。现已加 `gzip_ok()` 全量 CRC/ISIZE 校验 + "已完整则跳过"。

---

## 六、分步执行（每步独立可跑、有落盘产物、可中断）

时间预算：**机器合计 40–60 min**，全局峰值 < 1.5 GB。

| Step | 做什么 | 产物 | 机时 | 断点 |
|---|---|---|---|---|
| 0 | 环境自检（.Rlib/.venv/磁盘/内存/网络，11 项断言） | 体检表 | <1 min | — |
| 1 | 下载 CORE-8 + 参考系，校验尺寸/gzip/SHA-256 | `data/raw/`、`_manifest.tsv` | ~30 min（挂后台） | ✅ |
| 2 | 逐数据集读入 + 统一基因空间（SYMBOL 规范化、低表达过滤、剔除 `^mt-`）。**绝对禁止跨 GSE cbind/merge**（脚本内加 `stop()` 断言） | `data/interim/<GSE>_expr.rds` | ~5 min | ✅ |
| 3 | 逐数据集单独 DE（计数→limma-voom 或 DESeq2；芯片→limma）。固定列序 `gene, logFC, AveExpr, t, P.Value, adj.P.Val` | `results/de/<GSE>_de.csv` | ~10 min | ✅ |
| 4 | **第一幕** 跨模型 logFC 相似度矩阵 + 热图 + 聚类 | `results/act1/` | ~2 min | ✅ |
| 5 | **去卷积**：流式构建 S → NNLS 解 P（+ 80/20 基因级 held-out R²） | `data/ref/`、`results/deconvolution/` | ~12 min + 727 s 建 S | ✅ |
| 6 | **第二幕（核心）** MRM + 模型错配置换 | `results/act2/` | ~5 min | ✅ |
| 7 | **第三幕** 残差化 → 残差上重跑 DE 与跨模型比较（**残差必须用 limma 不能用 DESeq2，因为会出现负值**） | `results/act3/` | ~10 min | ✅ |
| 8 | 校准与敏感性（伪 bulk 已知比例对照） | `results/calibration/` | ~10 min | ✅ |
| 9 | 出图出表（fig01–fig15 + table01–table10） | `results/figures/`、`results/tables/` | ~5 min | ✅ |

**人时（真正的成本）**：元数据分组核对 1–2 h、看图写报告 4–8 h，合计 **8–15 h**。
**会话切分建议**：拆成 7 个会话 × 每周 2 小时 ≈ **3.5 周**（符合「初试前每周 2 小时、初试后到 3 月复试之间真正跑」的节奏）。

### Step 5 的验收判据（硬）
- `S` 必须是 **37,241 × 11**（行是基因）；P 每列和 == 1
- 重建 R² 中位数 > 0.6（**已在 GSE236388 上实测 0.680–0.742 ✅**）
- **ΔP 方向必须是 KO 组 germ-cell 下降、Somatic 上升** —— **这是整个项目的正向对照，已 PASS（−0.194 / +0.194）✅**
- **S 必须做 CPM 标准化**：实测未标准化的 S 会让**体细胞比例高估 0.26**，且偏倚方向恰好与待检验的组成位移同向。CPM 把误差从 0.151 降到 0.033（**4.6×**）。
  - S 的行和（raw mean counts/cell）：Meiotic 10,576 / Postmeiotic 6,058 / Somatic 2,418 / Spermatogonia 2,179 → Meiotic 是 Spermatogonia 的 **4.9 倍**，不标准化必然偏。

### 参考系的三个解析陷阱（已实测）
1. `DGE.txt.gz` 第 1 行**只有 34,633 个细胞条码、没有基因名列**；数据行有 34,634 字段。pandas 必须 `header=None, skiprows=1, index_col=0`，且**不能传标量 dtype**（会报 `could not convert string to float: '0610005C13Rik'`），要传 `{i: np.float32 for i in range(1, ncell+1)}`。
2. `PerCellAttributes` 是 **CRLF 行尾**、表头**自己以 `#` 开头**、文件**末尾无换行符** → `wc -l` 报 34636、`awk 'NR>5'` 报 34632，**都比真实 34633 少 1**。唯一可信计数是 `len(labels)`。
3. 该系列**混装全睾丸与富集/分选样本**（ST1–6 富集 none，其余 diploid/interstitial/sertoli/spermatogonia）→ **必须先用 PerCellAttributes 过滤**，否则参考系自身被分选污染，自相矛盾。

---

## 七、不确定性（必须在论文里如实写）

| 级别 | 问题 | 对策 |
|---|---|---|
| 🔴 U-1 | **检验力**：M=6–8 时若真实 compShare≈0.75，可能无法与 0.90 区分 | 用 MRM（不依赖比值分母）；报区间不报两位小数；Limitations 明写 |
| 🔴 U-2 | **参考系是成年小鼠**，而部分 bulk 是幼鼠（GSE38096 是 P0、GSE68284 是 P21）或不同品系 → **系统性偏倚而非随机噪声** | 报告里明写；按年龄分层做敏感性 |
| 🟠 U-3 | 节点标注只有 161/220 来自文献明确陈述；但只留可靠节点会进一步压小 M，与 U-1 直接冲突 | **需拍板取舍**（见下） |
| 🟡 U-5 | **校准有循环性**：用 GSE112393 合成伪 bulk 再用它去卷积必然成功，**不验证参考系匹配真实组织** | 必须补跨范式校准（如用分选参考 PMC4837615）；**预判很可能不通过** |
| 🟢 U-7 | ~~NNLS 正交性机械压低残差相似度~~ **方向是反的，已排除** | 真实风险是残差相似度被**抬高**（a=0 时同节点比跨节点高 **δ=+0.092, z=+2.57**）⇒ **第三幕必须把 δ 当假阳性零点减掉** |

### 已通过的 8 条校准判据（在看真实数据前固定）
C1 失配参考系比例回收 MAE < 0.05（实测 **0.024** ✅）/ C2 a=0 时 compShare ∈ [0.8,1.2]（实测 **1.020** ✅）/ C3 随真信号单调下降（1.02→0.82→0.50→0.21→0.07 ✅）/ C5 正交性不机械压低残差（差 **−0.005** ✅）/ C6 CPM 至少优于 raw 2×（实测 **4.6×** ✅）/ C7 置换零分布非退化 ✅ / C8 held-out R² 有判别力 ✅

---

## 八、需要你拍板的三件事

1. **U-3 vs U-1 的真实张力**：节点标注只有 161/220 来自文献明确陈述；但只留可靠节点会进一步压小 M —— 而 M=6–8 的检验力已经是最大短板。**这两条约束互相挤压，需要你决定取舍**（建议：主分析用可靠节点，全节点做敏感性）。
2. **Unknown 细胞（占 6.37%）不该默默丢弃**：实测丢弃比并入 Postmeiotic 的回收误差差 **2.4 倍**（0.0496 vs 0.0205），且偏倚方向与科学假设同向。**建议主分析丢弃 + 另跑"并入"臂作敏感性，不要默默丢。**
3. **Cal-1 的循环性**：用 GSE112393 合成再用它去卷积必然成功。**必须补跨范式校准**，且我预判很可能不通过 —— 如果不通过，结论强度要相应下调。

---

## 九、风险排序

| 级别 | 风险 | 处理 |
|---|---|---|
| 🔴 1 | **Step 5 稠密读入**（报错 `vector memory limit of 16.0 Gb reached` 看着像"内存不够"容易误判成"换机器"） | 用流式，峰值 47.7 MB |
| 🟠 2 | **元数据分组核对**（WT/KO 搞反**不会报错**但结论全反） | Step 2/3 加断言；人工逐条核对 1–2 h |
| 🟠 3 | `GSE244681_RAW.tar` **462 GB**（44 天） | 直接放弃 |
| 🟡 4 | 并行跑 8 个 DESeq2（8×2 GB = 16 GB 满） | 串行 |
| 🟡 5 | NCBI 503 限流 | 挂后台 + 块级续传 + 指数退避 |
| 🟡 6 | 矩阵朝向/对数空间搞混 | Step 2 加朝向断言（`grep` 基因名在行还是列） |
| 🟡 7 | GSE259417 只有作者 DEG | 已知：只能算 logFC，不进去卷积 |

---

## 十、文件地图

```
spermatogenic-arrest-plan/
├── README.md                    ← 仓库入口
├── PLAN.md                      ← 本文件（完整方案）
├── LICENSE
├── .gitignore
├── docs/
│   ├── 01-environment.md        (335 行) 环境实测：21 包逐包结果 + 内存峰值 + 做不到的事
│   ├── 02-datasets.md           (230 行) 数据集选择逻辑 + 库存纠错 + 下载预算
│   ├── 03-null-model.md         (843 行) 零模型规格 + 校准判据 + 不确定性清单
│   └── 04-pipeline-runbook.md   (913 行) 端到端 runbook（含每步失败处理）
├── code/
│   ├── 04_normalize_signature.py    S → CPM 标准化（含行和/列名逐字节校验）
│   ├── 05_deconvolution.R           NNLS 去卷积 + tissue_class 硬门槛 + held-out R²
│   └── 06_validate_deconvolution.py 端到端验证器 + 正向对照（已 PASS）
├── scripts/
│   ├── 01_geo_search.py … 03_curate.py   数据集盘点（生成 inventory.csv）
│   ├── build_signature.py                流式构建 S（峰值 47.7 MB）
│   └── fetch_geo.py                      并发 range 下载器（块级续传 + gzip 校验）
├── envs/                        setup.sh / install_pkgs.R / _drive.sh
├── data/
│   ├── README.md                原始数据的下载说明（数据本身不入库）
│   ├── inventory.csv            220 行自动盘点（⚠️ 有 15 处标错，见 §四）
│   ├── dataset_selection.csv    57 行 × 19 列（core 8 / extended 4 / rejected 30 / mislabel 15）
│   └── samples.tsv              86 行逐样本（含 tissue_class，去卷积脚本读它）
└── notes/
    ├── 原计划书.md              ← 参考（已不采用，见 §二）
    ├── 思路整理.md              ← 科学问题来源（已采用）
    └── 数据集清单.md
```

运行时生成、不入库的路径：`.Rlib/`（171 包 / 934 MB）、`.venv/`、`.cache/`、
`logs/`、`data/raw/`（原始表达数据）、`data/ref/`、`data/interim/`、
`envs/vendor/`、`envs/probe/`。

---

## 十一、马上可以做的第一条命令

```bash
cd spermatogenic-arrest-plan

# 0) 先把两个数据文件放好（原始数据不入库，见 data/README.md）
#    data/ref/signature_GSE112393_cpm.csv.gz
#    data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz

# 1) 环境自检（幂等，约 15 s）
bash envs/setup.sh

# 2) 复现已经跑通的正向对照（证明地基可用，约 1 min）
MPLCONFIGDIR="$PWD/.mplcache" .venv/bin/python code/06_validate_deconvolution.py \
    --sig data/ref/signature_GSE112393_cpm.csv.gz \
    --bulk data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz \
    --wt-pattern WT --ko-pattern cKO
# 期望输出：POSITIVE CONTROL PASS，germ-cell delta = -0.194，R² 0.680–0.742
```
