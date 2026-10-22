# 数据集与去卷积参考系：最小可行数据集清单

> 历史盘点记录，非当前方法或纳入状态的权威说明。下文将 GSE38096 称为双色芯片的记录有误：实际为单 Cy3 通道。Hnrnpgt 旧转录本映射、TPM 与 counts 的处理、样本数量及参考可用性以 [当前方法](METHODS.md)、[当前结果](RESULTS.md) 和 `data/analysis_assets.json` 为准。旧的“能去卷积”不等于比例经过验证。

> 产出人：dataset-curator（task-2）　日期：2026-10-16
> 配套文件：`data/dataset_selection.csv`（57 行，入选/淘汰/库存纠错）、`data/samples.tsv`（85 行，逐样本清单）

---

## 0. 本报告最重要的一条（老师最可能追问的点）

**去卷积只能用于"全睾丸（whole testis）"样本。**

分选/富集细胞群体的**组成是实验者用分选决定的，不是测量出来的**；对这类样本做去卷积，等于用参考系去"推断"一个你早就知道答案、而且答案由操作定义的量——结果没有意义，且会自我循环。

因此本篇对**每一个**数据集都先判定 `tissue_class`，只有 `whole_testis` 才允许进入去卷积。判定依据不是 `data/inventory.csv` 的 `tissue` 字段（该字段 220 行里有 53 行是"待查"，且我们发现多处标错），而是 GEO 的 **`targ=gsm` 逐样本记录**（`!Sample_source_name_ch1` + `!Sample_characteristics_ch1` 的 `tissue:` 行）与 `!Series_overall_design` 原文。所有判定都附了可复查的原文引用。

`tissue_class` 三值：`whole_testis` / `sorted`（分选、富集、分级、免疫沉淀等）/ `scrna`（单细胞）。

---

## 1. 入选清单（CORE-8，覆盖三个节点）

| # | GSE | 模型 | 节点 | tissue_class | n (KO/WT) | 年龄 | 平台 | supp 类型 | 可去卷积 |
|---|-----|------|------|--------------|-----------|------|------|-----------|----------|
| 1 | **GSE38096** | Sin3a 生殖细胞特异 cKO (Ddx4/Vasa-cre) | 精原细胞/起始 | whole_testis | 6 (2/2 + 2 HET) | P0 | GPL13912 Agilent 4x44K | RAW.tar 每样本 txt | yes（需归一化） |
| 2 | **GSE42601** | Dazap1-Fn 突变（减数分裂 I 前阻滞） | 精原细胞/起始 | whole_testis | 6 (3/3) | P21 | GPL1261 Mouse430_2 | RAW.tar 每样本 CEL | yes（需 RMA） |
| 3 | **GSE294209** | Rnu11 cKO（次要剪接体 U11 snRNA） | 减数分裂前期 | whole_testis | 6 (3/3) | P14 | GPL24247 NovaSeq 6000 | RAW.tar 每样本 gene/iso TPM | yes（TPM，非 counts） |
| 4 | **GSE236388** | XAP5 cKO（**仅 D16 睾丸臂**） | 减数分裂前期 | whole_testis | 6 (3/3) | D16 | GPL24247 NovaSeq 6000 | counts 矩阵 | **yes（首选，原始 counts）** |
| 5 | **GSE101511** | Hnrnpgt Cre-Lox cKO | 减数分裂前期 | whole_testis | 6 (3/3) | 18 dpp | GPL19057 NextSeq 500 | RAW.tar 每样本 counts | yes（需转录本→基因汇总） |
| 6 | **GSE283904** | N6amt1-SKO（生殖细胞特异） | 减数分裂前期 | whole_testis | 8 (4/4) | 成年 | GPL24247 NovaSeq 6000 | counts 矩阵 | **yes（首选，原始 counts）** |
| 7 | **GSE68283** | Rfx2-/- （圆形精子期末完全阻滞） | 减数分裂后/精子形成 | whole_testis | 9 (4/5) | P29-30 | GPL17021 HiSeq 2500 | raw counts 矩阵 | **yes（首选，原始 counts）** |
| 8 | **GSE68284** | Rfx2-/- （同论文早期时间点） | 减数分裂后/精子形成 | whole_testis | 11 (6/5) | P21 | GPL17021 HiSeq 2500 | raw counts 矩阵 | **yes（首选，原始 counts）** |

节点覆盖：**精原细胞/起始 2 个**（GSE38096、GSE42601）、**减数分裂前期 4 个**（GSE294209、GSE236388、GSE101511、GSE283904）、**减数分裂后/精子形成 2 个**（GSE68283、GSE68284）。三节点全覆盖。

### 每个数据集的"全睾丸"证据（原文）

- **GSE38096** — GSM `source_name_ch1` 字面为 `whole testes, Vasa-cre;Sin3a<deleted/floxed> mice (VSKO)`；`characteristics_ch1 = tissue: whole testis`（6/6 样本）；`overall_design`: *"Gene expression in the **whole testes** of age-matched Ddx4(Vasa/Mvh)-cre mediated conditional Sin3a gene-targeted(VSKO) mice, germ cell-specific Sin3a hemizygous (HEMI) mice, and wild type (WT) mice was measured at postnatal day 0."*
  注意：系列为 **VSKO / HEMI / WT 三组各 2 只**，只有 VSKO-vs-WT 是干净的 2v2 KO/WT 对比；P0 为新生期睾丸（尚无减数分裂后细胞），参考系必须按年龄匹配。
- **GSE42601** — GSM `source_name_ch1` = `testes, wild-type` / `testes, Dazap1 mutant`，`tissue: testes`，`age: P21`；`overall_design`: *"We used wild-type and Dazap1 mutant mouse testes. Each genotype has three replicates. 6 total samples were analyzed."*
- **GSE294209** — GSM `characteristics_ch1 = **tissue: Whole testis**`（6/6）；`overall_design`: *"RNA-seq profiling of postnatal day 14 **whole-testes** dissected from control and mutant (Rnu11-cKO) animals (n=3)."*
- **GSE236388** — 可用的 6 个 GSM 标题为 `Testes, WT, D16, rep1-3`（GSM7529027-29）与 `Testes, XAP5 cKO, D16, rep1-3`（GSM7529030-32），`tissue: testes`；`overall_design` 明确 *"performed RNA-seq analysis using **D16 WT and XAP5 cKO testes** samples. Three biological replicates were included for each genotype."*
- **GSE101511** — GSM `source_name_ch1 = Testes`，`tissue: Testes`（6/6）；`overall_design`: *"Gene expression profiling by RNA-Seq of mouse **testes** 18 days post-partum. Samples from C57BL/6 background, either wild type (n=3) or HnRNPGT Cre-Lox knockout (n=3)."* 基因型：Hnrnpgt-/- = 293/320/263（GSM2705201/02/04）；wild type = 321/279/288（GSM2705203/05/06）。
- **GSE283904** — GSM `source_name_ch1 = testis`，`tissue: testis`（8/8）：CTL-1..4 = GSM8674146-49，N6A-SKO-1..4 = GSM8674150-53；`overall_design`: *"...**testes** from the CTL and N6amt1-SKO mice were sent to RNA-seq analysis."*
- **GSE68283 / GSE68284** — GSM `source_name_ch1 = **Decapsulated testis**`，`tissue: Testis`，`age: post natal day 29-30`（GSE68283）/ P21（GSE68284），`genotype: wild type` 或 `Rfx2-/-`。

### supp 文件结构与"是否真能去卷积"（已逐文件验证）

| GSE | 文件 | 字节 | 表头/结构 | 判定 |
|-----|------|------|-----------|------|
| GSE236388 | `GSE236388_XAP5_counts.txt.gz` | 285,621 | `GeneName ⇥ XAP5 WT replicate1-3 ⇥ XAP5 cKO replicate1-3`，23,502 基因 | ✅ 原始 counts |
| GSE283904 | `GSE283904_Gene_count.txt.gz` | 347,406 | `gene_id ⇥ N6_CTL_1-4 ⇥ N6_KO_1-4`，17,121 基因 | ✅ 原始 counts |
| GSE68283 | `GSE68283_Raw_counts.txt.gz` | 389,507 | `gene_id ⇥ 141KO 173WT 174WT 182KO 183KO 194WT 197WT 207WT 209KO`，23,074 基因 | ✅ 原始 counts |
| GSE68284 | `GSE68284_Raw_counts.txt.gz` | 474,203 | `gene_id ⇥ KO10 KO11 KO6 KO7 KO8 KO9 WT1..WT5`，23,074 基因 | ✅ 原始 counts |
| GSE294209 | `GSE294209_RAW.tar` | 10,700,800 | 12 个成员 = GSM8900412-17 × {gene_tpm, iso_tpm}_estimates.txt.gz | ⚠️ **每样本 TPM**（可重建矩阵，但非 counts） |
| GSE101511 | `GSE101511_RAW.tar` | 1,914,880 | 6 个 `GSM*_*.count.txt.gz`，ENSMUST 转录本级 2 列 | ⚠️ 需先汇总到基因水平 |
| GSE235777(备) | `GSE235777_RAW.tar` | 2,027,520 | 6 个 `GSM7509320-25_*ReadsPerGene.out.tab.gz`（STAR） | ⚠️ 取第 2 列 |
| GSE38096 | `GSE38096_RAW.tar` | 83,292,160 | 每样本 Agilent FE txt.gz（双色芯片） | ⚠️ 需背景校正 + loess |
| GSE42601 | `GSE42601_RAW.tar` | 24,145,920 | 6 个 `GSM1045752-57_*.CEL.gz`（GSM1045752_1312W.CEL.gz … GSM1045757_1364F.CEL.gz） | ⚠️ 需 RMA |
| GSE213442(备) | `GSE213442_cpm_matrix.txt.gz` | 1,415,035 | `control1-4 kpn1-4`，24,988 基因 | ⚠️ CPM，非 counts |
| GSE214201(备) | `GSE214201_cpm_matrix.txt.gz` | 975,523 | `control1-4 iqch1 iqch2 iqch4`，17,128 基因 | ❌ **只有 3 个 KO 列**（缺 iqch3） |

---

## 2. GSE259417：只能算 logFC，**不能去卷积**（硬结论）

GEO 上该系列的 supp **只有** `GSE259417_P16_DEG.csv.gz`（2,400,869 B，与 FTP 列出的 2.3M 一致）。结构已逐列验证：

```
['', 'Row.names', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj',
 'WT_rep1', 'WT_rep2', 'KO_rep1', 'KO_rep3', 'chrom', 'type', 'symbol']   29,332 数据行
```

这是作者用 DESeq2 算好的 DEG 表：既有 logFC/padj，也带了 **2/3 个 WT 和 2/3 个 KO** 的归一化计数列（注意 `KO_rep2` 缺失，命名也不规则）。

**结论：能算 logFC ✅，不能做去卷积 ❌。** 去卷积需要**每个样本**完整基因表达向量以估计组成比例；这里既没有全基因组 × 全样本矩阵，样本数也只有 2v2。该数据集只能作为 logFC 层面的外部一致性验证集，不能进入 Act2 的组成零模型。

> 库存纠错：`data/inventory.csv` 中 GSE259417 的 `supp_state = matrix` **不正确**，应记为 `author_DEG_only`。

---

## 3. GSE112393 参考系核验结论（去卷积参考）

**结论：可用，且方向与任务描述完全一致，但必须按 `PerCellAttributes` 过滤，不能整份直接当组成参考。**

系列：Green et al., *A Comprehensive Roadmap of Murine Spermatogenesis Defined by Single-Cell RNA-Seq*, Dev Cell 2018；PMID **30146481**。tissue_class = **scrna**（单细胞，34,633 细胞），**永远不作为 bulk 样本**，只作为去卷积的**特征矩阵（signature matrix）**来源。

已下载并核验：

| 文件 | 字节 | 核验结果 |
|------|------|----------|
| `GSE112393_MergedAdultMouseST25_DGE.txt.gz` | 81,261,540 | `gzip -t` 通过（完整） |
| `GSE112393_MergedAdultMouseST25_PerCellAttributes.txt.gz` | 1,495,253 | 已下载 |

**矩阵方向（已实测，非推测）**：第 1 行 = **34,633 个细胞条码**（如 `ST1_ATCTTGCACATC`、`ST1_CCGCCTTAGCTN`…），**没有基因名列**；从第 2 行起第 1 列才是基因名（第 2 行以 `0610005C13Rik  0  0` 开头），该行共 34,634 个字段 = 1 基因名 + 34,633 计数。即 **行 = 基因、列 = 细胞**，与任务描述一致（"第1行是细胞条码 34633 列，第1列是基因名"）。

**`PerCellAttributes.txt.gz` 可用性**：列 1 = CellBarcode、列 2 = Dataset（样本标签，如 ST1…ST25）、列 3 = CellType。**列 3 是 1–11 的数值编码，不是字符串**，下游脚本必须先查表再比较，否则会把 `"9"` 当字符串而匹配失败：

```
1-InnateLymph  2-Macrophage  3-Endothelial  4-Myoid  5-Leydig  6-Sertoli
7-Unknown      8-SPG         9-Scytes       10-STids   11-Elongating
```

**关键警告**：GSE112393 的 25 个样本里，**同时包含全睾丸样本（ST1–ST6，`cell type enrichment: none`）和经过富集/分选的样本**（如 `diploid`、`interstitial`、`sertoli`、`spermatogonia`，以及"去除了单倍体细胞"的富集）。若只想要"全睾丸来源的细胞组成参考"，**必须先用 `PerCellAttributes` 的 Dataset / enrichment 标签把富集样本过滤掉**，否则参考系本身就被分选"污染"，与第 0 节的原则自相矛盾。

**回退方案**：PMC4837615 的 FACS 纯化阶段特异生殖细胞可作为**特征基因签名（signature genes）**的备选来源；但它本身是分选数据，**只能用来定义 marker 基因列表，不能用来定义组成比例**。Act2 的组成零模型必须以 GSE112393 的 whole-testis 子集为准。GSE120506（人睾丸图谱）为人类对照，不参与小鼠主线。

---

## 4. 被淘汰者及理由

### 4.1 库存标成"全睾丸"，实际不是（**重点纠错**）

| GSE | 库存标签 | 实际情况（GEO 原文） |
|-----|----------|---------------------|
| **GSE323391** | tissue=`全睾丸` | 系列标题即为 *"...in mouse testicular **organoids**"* → **体外类器官**，非体内全睾丸 |
| **GSE237030** | tissue=`全睾丸` | *"in vivo gonadal cells and ES cell-derived gonadal somatic cells were analyzed using **10X Chromium Single Cell**"* → 单细胞 + ES 诱导体细胞 |
| **GSE125495** | tissue=`全睾丸` | GSM `source_name_ch1 = **mouse female gonadal soma**` ×10 → **雌性**性腺体细胞，非睾丸 |
| **GSE310340** | 归入激素/内分泌（暗示睾丸相关） | GSM `source_name_ch1 = **Kidney** ×10 + **Liver** ×10` → 与睾丸无关 |
| **GSE23119** | tissue=`全睾丸` | *"isolated **pure populations of spermatogonia**…"*；提取方案原文：*"germ cells, in suspension, were separated by sedimentation with use of a 2-4% BSA gradient … Purity of spermatogonia was estimated and was routinely higher than **90%**."*（STATPUT 分选）→ **分选精原细胞** |
| **GSE158471** | tissue=`待查` | GSM `source_name_ch1 = LZ (leptotene and zygotene)` ×4 / `PD (pachytene and diplotene)` ×4，`cell type: Spermatocytes` → **分级分选精母细胞** |
| **GSE253154** | tissue=`待查` | `overall_design`: *"...we **sorted** differentiated spermatogonia and conducted RNA-seq."* → 分选 |
| **GSE288152** | tissue=`待查` | `overall_design`: *"Spermatogenic cells were **sorted** from 8-week-old wildtype and Tex10f/-; Stra8iCre mice."* → 分选 |
| **GSE324232** | tissue=`待查` | `overall_design`: *"bulk RNA-seq for (1) … pachytene spermatocytes (PS); (2) … round spermatids (RS)"* → **分选 PS/RS** |
| **GSE236388** | tissue=`全睾丸`, n=12 | 系列**混装**：6 个 D16 睾丸样本（全睾丸 ✅）+ 6 个 D56 附睾尾**精子**样本（`tissue: sperm`，`source_name_ch1 = sperm`）❌。库存 n=12 把两者合并了 |
| **GSE221433** | tissue=`全睾丸` | 全睾丸判定正确，但模型是**长期西方饮食**（18 周）而非遗传性阻滞 → 模型类别不同，不能当"遗传阻滞共性"的证据 |

### 4.2 测定类型不是表达矩阵（无基因水平表达）

GSE227944（FIGNL1/RAD51 ChIP-SSDS）、GSE113710（Dazl HITS-CLIP）、GSE181921（TH2A/TH2B 占据）、GSE206536（SRSF2 LACE-seq）、GSE53915（RNF17 5'RACE）、GSE57241（RNF17 RIP-seq）、GSE67683（MIWI/piRNA 非编码）、GSE68282（GSE63130 的 ChIP 臂）、GSE197610（免疫沉淀 RNA）、GSE139399（**组织确为全睾丸**，但 supp 只有核糖体印迹的 translation-efficiency 表，无 counts）。

### 4.3 分选 / 单细胞 / 混合模态

GSE230747（preleptotene 富集，SuperSeries）、GSE244681（Hi-C+ChIP+RNA 混合，RAW.tar **462 GB**）、GSE228454（sorted + scRNA，RAW.tar 1.9 GB）、GSE70889（分选生殖细胞）、GSE55180（prospermatogonia 分选）、GSE290629（6 个分选阶段）、GSE252121（Thy1+/c-Kit+ 分选）、GSE109033（10X）、GSE222626（圆形精子）、GSE318472（精原细胞）、GSE274959（scRNA + 基因型 multiplexed，RAW.tar 743 MB）、GSE95138（分选精子/精母细胞）、GSE174731（scRNA，n=2）。

### 4.4 非睾丸组织 / 体外

GSE238226（卵母细胞 + 胚胎）、GSE273732（E14Tg2a 细胞系）、GSE82002（GV 卵母细胞）。

### 4.5 全睾丸但被降级为 EXTENDED / 备选

- **GSE213442 (Kpna2)**、**GSE214201 (Iqch)**：全睾丸干净（`source_name_ch1 = Adult testis`），CPM 矩阵；但 GSE214201 矩阵**只有 3 个 KO 列**（`iqch1 iqch2 iqch4`，缺 iqch3），与系列登记的 4 个 KO 不符。
- **GSE235777 (Msl3)**：全睾丸 P45，3 KO / 3 WT（**修正库存的 3/6 错误**），counts 可重建；但论文标题即 *"Loss of function of male-specific lethal 3 (Msl3) **does not affect spermatogenesis** in rodents"* → **表型阴性模型**，只能当阴性对照，不能当阻滞案例。
- **GSE259417 (Hsf5)**：全睾丸 P16，但只有作者 DEG → 见第 2 节，logFC-only。
- **GSE57745 (Dnmt3L)**：全睾丸，但 supp 内容**不是基因计数矩阵**——6 个 `*_RNA.txt.gz` 的行是转座子/重复家族（LTR64、Charlie2a、L1MB7…），仅 781 行；且每基因型每时间点 n=1，两个时间点还用 HET 作对照。
- **GSE114042 (Dazl)**：全睾丸但 xlsx 是作者 DEG 表（仅 log2FoldChange/Pvalue/FDR），且 n=1 Wt vs 1 SKO。
- **GSE121326 (RALDH/ATRA)**：全睾丸已确认，但论文结论是**阴性**（*"the GC-derived source of ATRA does not perform any specific role"`）→ 无表型。
- **GSE3436 (Spo11-/-)**：全睾丸已确认（`15d Spo11-/- whole testis` ×4、`12d wild type whole testis` ×3、`12d Spo11-/- whole testis` ×2），但 2006 年 **spotted cDNA 芯片**，无现成 counts → 只能作遗留 logFC 验证。

---

## 5. 待办（已解决）

### 5.1 GSE3436 平台 = **GPL977，NIA15K，technology = spotted DNA/cDNA**
**不是** Affymetrix。因此只能算 logFC/比值，**去卷积证据弱**（需先做背景校正与比值归一化，再映射到基因水平）。

### 5.2 GSE23119 平台 = **GPL1261，[Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array**，technology = in situ oligonucleotide
（注：另有同名 GPL23119 = SOMAscan 蛋白检测，属同名巧合，与本系列无关。）supp 仅 RAW.tar + filelist.txt（CEL 级别）→ 需 RMA，**同属"能算 logFC 但去卷积证据弱"档；且该数据集本身是分选精原细胞，已被否决。**

### 5.3 PMID 31037746（2019 Andrologia）——**全文未取得**

- 元数据（NCBI E-utilities esummary 实测）：**Kui F, Ye H, Chen XL, Zhang J. "Microarray meta-analysis identifies candidate genes for human spermatogenic arrest." *Andrologia* 51(8):e13301, 2019 Sep.** DOI **10.1111/and.13301**。
- 可及性：PubMed 的 `articleids` 中 **PMC 列表为空**；Europe PMC REST 返回 `isOpenAccess: N, inEPMC: N, pmcid: None` → **无任何开放获取全文副本**，Wiley 付费墙。
- **因此：本文全文未取得，本报告只依据其摘要，不臆测其正文内容。** 摘要要点：178 例人睾丸样本、7 个芯片数据集、按 Johnsen 评分分 4 种病理表型 + 正常，分表型做 DEG meta 分析；"G1 期有丝分裂细胞周期阻滞是减数分裂前阻滞的显著特征"；15 个 p53 互作蛋白（含 ABL1、HDAC2）。
- **对"该文是否讨论细胞组成"这一问题的回答：摘要中没有任何关于细胞类型组成/去卷积的讨论**；正文无法获取，故**无法确认**其是否在正文中处理了组成混杂——这一空白应作为我们相对该 precedents 的增量点，但不得声称该文"忽略了组成"。

---

## 6. 下载体积与折算小时数

按任务口径 **126 KB/s** 折算（同时给出实测 105 KB/s 的值作对照——本机 8 并发 range 实测约 98–107 KB/s）。

| 口径 | 字节 | MB | MiB | 折合小时 @126 KB/s | @105 KB/s（实测） |
|------|------|----|-----|--------------------|-------------------|
| **CORE-8 入选 supp** | 121,550,497 | 121.6 | 115.9 | **0.26 h**（≈16 分钟） | 0.31 h |
| **GSE112393 参考系**（DGE + PerCellAttributes） | 82,756,793 | 82.8 | 78.9 | **0.18 h**（≈11 分钟） | 0.21 h |
| **CORE-8 + 参考系 合计** | **204,307,290** | **204.3** | **194.8** | **0.44 h**（≈27 分钟） | 0.53 h |

明细（CORE-8）：GSE38096 83,292,160 · GSE42601 24,145,920 · GSE294209 10,700,800 · GSE101511 1,914,880 · GSE68284 474,203 · GSE68283 389,507 · GSE283904 347,406 · GSE236388 285,621。

参考系明细：GSE112393_DGE.txt.gz 81,261,540 · GSE112393_PerCellAttributes.txt.gz 1,495,253。

**一句话**：全部必需数据 (<205 MB) 在 126 KB/s 下约需 **27 分钟**；即使按实测 105 KB/s 也只有约 **32 分钟**。这与最初担心的"数百 GB"（如 GSE244681 的 462 GB）相差三个数量级——**下载预算从来不是本项目的瓶颈，唯一的瓶颈是找到足够多"确为全睾丸"的阻滞模型**。

> 下载状态说明：GSE38096_RAW.tar 目前磁盘上为 **24,543,601 B（部分文件，完整应为 83,292,160 B）**，其为 8 个核心数据集中最后完成的一个；其余 7 个已完整落盘并通过 `gzip -t`/`tar -tf` 校验。上表按**完整体积**计。

---

## 7. 库存纠错清单（待定是否回写 inventory.csv）

| GSE | 库存字段 | 正确值 |
|-----|----------|--------|
| GSE323391 | tissue=全睾丸 | 睾丸**类器官**（体外） |
| GSE237030 | tissue=全睾丸 | scRNA-seq（ES 诱导 + 体内体细胞） |
| GSE125495 | tissue=全睾丸 | **雌性**性腺体细胞 |
| GSE310340 | 归入睾丸相关 | Kidney + Liver |
| GSE259417 | supp_state=matrix | 仅作者 DESeq2 DEG，**无 counts** |
| GSE114042 | supp_state=matrix | 仅作者 DEG 表，**无 counts** |
| GSE57745 | supp_state=matrix | supp 为 781 行转座子家族计数（非基因水平） |
| GSE236388 | tissue=全睾丸, n=12 | 仅 6 个为全睾丸（D16）；另 6 个为附睾尾精子（D56） |
| GSE235777 | n_ko/n_wt = 3/6 | 实为 **3 KO / 3 WT** |
| GSE101511 | n_ko/n_wt 为空 (0/0) | 实为 **3 KO / 3 WT** |
| GSE158471 | tissue=待查 | 分选 LZ/PD 精母细胞 |
| GSE253154 | tissue=待查 | 分选精原细胞 |
| GSE288152 | tissue=待查 | 分选生精细胞 |
| GSE324232 | tissue=待查 | 分选 PS + RS |
| GSE221433 | tissue=全睾丸（正确） | 但模型为**西方饮食**，非遗传阻滞 |

---

## 8. 给下游（Act1/Act2/Act3）的接口约定

- `data/dataset_selection.csv` 的 **`tissue_class`** 列是去卷积的**硬门槛**：`code/05_deconvolution.R` 应直接拒绝任何 `tissue_class != "whole_testis"` 的输入。
- `data/dataset_selection.csv` 的 **`deconvolvable`** 列区分 `yes` / `yes (with caveat)` / `partial` / `no`——"yes with caveat" 表示需要先做归一化/汇总/RMA 等预处理才能进入去卷积。
- `data/samples.tsv`（85 行，列：`gse, sample_id, title, group, tissue_class, node, age, platform, note`）是逐样本的 DE 与去卷积输入清单。其中：
  - `tissue_class = whole_testis` 共 **79** 个样本 → 可用于 DE 与去卷积；
  - `tissue_class = sperm_EXCLUDED` 共 **6** 个样本（GSE236388 的 D56 精子臂）→ **保留在文件中但显式标记排除**，防止被误当成全睾丸并入。
- 尚不确定项（本报告内已显式标注，未留空白）：GSE294209 与 GSE283904 的 `!Series_pubmed_id` 在 GEO 记录中缺失，其 PMID 需另行从论文核对（标记为 `NEEDS_REVIEW`）；GSE38096 的完整 tar 尚在下载中。

---

## 9. 引用（均为真实可查）

GEO 系列：GSE38096 · GSE42601 · GSE294209 · GSE236388 · GSE101511 · GSE283904 · GSE68283 · GSE68284 · GSE213442 · GSE214201 · GSE235777 · GSE259417 · GSE23119 · GSE158471 · GSE253154 · GSE288152 · GSE324232 · GSE323391 · GSE237030 · GSE125495 · GSE112393 · GSE3436 · GSE57745 · GSE114042 · GSE121326。

平台：GPL13912 · GPL1261 · GPL24247 · GPL19057 · GPL17021 · GPL977 · GPL18635。

文献 PMID（均经 NCBI E-utilities 核验）：
- **30146481** — Green CD et al., *A Comprehensive Roadmap of Murine Spermatogenesis Defined by Single-Cell RNA-Seq*, Dev Cell 2018（GSE112393 参考系）
- **31037746** — Kui F et al., *Microarray meta-analysis identifies candidate genes for human spermatogenic arrest*, Andrologia 2019;51(8):e13301, DOI 10.1111/and.13301（**全文未取得**）
- **23085237**（GSE38096）· **23965306**（GSE42601）· **30674417**（GSE101511）· **26162102**（GSE68283/68284）· **37456838**（GSE213442）· **36975404**（GSE214201）· **37847071**（GSE235777）· **38007088**（GSE221433）· **40009443**（GSE236388）· **23253600**（GSE23119）· **27149259**（GSE55180）
