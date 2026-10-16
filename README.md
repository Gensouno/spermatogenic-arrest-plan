# 生精阻滞转录组共性分析

小鼠生精阻滞（spermatogenic arrest）模型的全睾丸 bulk 转录组共性分析。核心问题是：

> 不同阻滞模型之间那套"共同差异基因"，有多少是真正共享的调控信号，
> 有多少只是**细胞组成变化**的副产物？

本仓库是完整的分析方案、方法学规格与可复现代码。

---

## 背景

睾丸 bulk RNA-seq 里被反复报告的"共识信号"（`Prm1/Prm2/Tnp1/Tnp2` 下调、
`Sox9/Wt1/Gata4/Amh/Cyp11a1` 上调）在生物学上就是「丢失了晚期生精细胞」这一件事的
另一种说法。任何造成生精阻滞的模型都会让晚期细胞减少、体细胞相对富集，
因此这类"共性"在多大程度上只是组成镜像，需要被定量检验，而不是被默认接受。

本方案采纳的问题框架是三幕结构：

| 幕 | 做什么 |
|---|---|
| 第一幕 | 天真版：跨模型 logFC 收敛，复现文献常见观察 |
| **第二幕（核心）** | 构造**只含细胞组成**的零模型，量化组成能解释多少跨模型相似度 |
| 第三幕 | 扣掉组成预测后，在残差里找真正剩下的信号 |

---

## 已经验证通过的地基

方案不是纸上推演。用真实数据做端到端去卷积，**项目自己的正向对照已经 PASS**：

```
数据集：GSE236388（XAP5 cKO 睾丸臂，D16，3 KO / 3 WT）
参考系：GSE112393（成年小鼠睾丸，37,241 genes × 34,633 cells → 4 类细胞）
```

| 指标 | 结果 |
|---|---|
| 重建 R² | 0.680 – 0.742（阈值 0.6） |
| germ-cell 比例 ΔP | **−0.194**（KO 下降） |
| Somatic 比例 ΔP | **+0.194**（KO 上升） |

一个生精阻滞敲除应当表现为 germ-cell 比例下降、体细胞比例相应上升——
方向与幅度都对。复现命令见文末。

---

## 快速开始

```bash
git clone <repo-url>
cd spermatogenic-arrest-plan

# 1) 环境自检（幂等，约 15 s）
bash envs/setup.sh

# 2) 复现已跑通的正向对照（约 1 min）
#    需要先下载两个数据文件，见 data/README.md
MPLCONFIGDIR="$PWD/.mplcache" .venv/bin/python code/06_validate_deconvolution.py \
    --sig  data/ref/signature_GSE112393_cpm.csv.gz \
    --bulk data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz \
    --wt-pattern WT --ko-pattern cKO
```

环境要求：macOS aarch64（Apple Silicon）、R ≥ 4.4、Python ≥ 3.10、**16 GB 内存**即可。
`envs/setup.sh` 会把 R 包装进项目内 `.Rlib`、Python 包装进 `.venv`，
**不写系统目录**。若你的 R 全局库可写，可自行简化为系统安装。

> GEO 下载受 NCBI 限速（单流 15–35 KB/s），仓库提供 8 并发分块下载器
> `scripts/fetch_geo.py`（实测约 6 倍加速、支持断点续传）。
> 原始数据不入库，见 [`data/README.md`](data/README.md)。

---

## 目录结构

```
.
├── PLAN.md                  总方案（建议从这里读起）
├── docs/
│   ├── 01-environment.md    本机环境实测：逐包安装结果、内存峰值、已知限制
│   ├── 02-datasets.md       数据集筛选：CORE-8 全套证据链 + 公开数据纠错
│   ├── 03-null-model.md     零模型的数学规格、校准判据、不确定性清单
│   └── 04-pipeline-runbook.md  端到端执行规格（Step 0–9，含每步失败处理）
├── code/
│   ├── 04_normalize_signature.py    参考谱 CPM 标准化
│   ├── 05_deconvolution.R           NNLS 去卷积 + 硬门槛 + held-out R²
│   └── 06_validate_deconvolution.py 端到端验证器 + 正向对照
├── scripts/
│   ├── 01_geo_search.py … 03_curate.py  数据集盘点
│   ├── build_signature.py               流式构建细胞类型特征谱
│   └── fetch_geo.py                     并发分块下载器
├── envs/setup.sh            幂等环境搭建
├── data/                    数据集清单与逐样本分组表（原始表达数据不入库）
└── notes/                   本项目的两版早期构想与数据集盘点全文
```

---

## 数据

分析使用 **8 个全睾丸** 生精阻滞模型（每个都逐样本核验过 GEO 原文）：

| GSE | 模型 | 阻滞节点 | n (KO/WT) | 年龄 | 平台 |
|---|---|---|---|---|---|
| GSE38096 | Sin3a 生殖细胞特异 cKO | 精原细胞/起始 | 6 (2/2+2HET) | P0 | GPL13912 |
| GSE42601 | Dazap1 突变 | 精原细胞/起始 | 6 (3/3) | P21 | GPL1261 |
| GSE294209 | Rnu11 cKO | 减数分裂前期 | 6 (3/3) | P14 | GPL24247 |
| GSE236388 | XAP5 cKO（仅 D16 睾丸臂） | 减数分裂前期 | 6 (3/3) | D16 | GPL24247 |
| GSE101511 | Hnrnpgt cKO | 减数分裂前期 | 6 (3/3) | 18 dpp | GPL19057 |
| GSE283904 | N6amt1-SKO | 减数分裂前期 | 8 (4/4) | 成年 | GPL24247 |
| GSE68283 | Rfx2 −/− | 减数分裂后/精子形成 | 9 (4/5) | P29–30 | GPL17021 |
| GSE68284 | Rfx2 −/−（同论文 P21） | 减数分裂后/精子形成 | 11 (6/5) | P21 | GPL17021 |

去卷积参考系：**GSE112393**（成年小鼠睾丸单细胞图谱，Green et al., *Dev Cell* 2018）。

> **去卷积只适用于全睾丸样本。** 分选细胞群体的组成由分选决定、不是测量出来的，
> 其比例不是生物学估计。`code/05_deconvolution.R` 对此有硬性 `stop()` 守卫。

⚠️ 公开数据集标签有相当多的错误。`data/inventory.csv` 是自动盘点结果，
其中 **15 条经逐条核验为误标**（睾丸类器官被当作全睾丸、附睾尾精子混入睾丸臂、
n 值记错等）。**请以 `data/dataset_selection.csv` 为准**，
纠错明细见 [`docs/02-datasets.md`](docs/02-datasets.md)。

---

## 方法要点

零模型预测的 logFC —— 组成位移经参考谱推进到表达空间，且必须带 log：

$$\hat c_m(g)=\log\!\big((S\,\bar p^{KO}_m)(g)+1\big)-\log\!\big((S\,\bar p^{WT}_m)(g)+1\big)$$

含义：**若每个基因的表达本身毫无变化、只有细胞组成变了**，模型 *m* 会「看起来」
有哪些基因上调/下调。三个用模拟实测出来的坑（详见
[`docs/03-null-model.md`](docs/03-null-model.md)）：

1. **最自然的估计量是错的。** 直接用「logFC 投影到 *S* 的列空间」作组成预测，
   *即使给完美参考系也只能报 0.585（真值 1.0）*，系统性低估约 40%——
   因为 log 是凹函数，组成效应不落在 *S* 的列空间里。
2. **朴素 Mantel 检验无信息。** 在「纯伪影」与「有真信号」两种相反情形下给出
   同一个 *p* = 0.158。改用**以模型对为单位的 MRM 回归**。
3. **置换单位必须是「模型」**，且必须用模型错配。基因打乱法的零分布
   sd = 0.007（过窄、大量假阳性）；模型错配法 sd = 0.129 才是正确的。

**检验力的诚实边界**：M = 6–8 时能回答「组成是不是主要解释」，
**不能**回答「37% 还是 44%」。因此 compShare 一律报点估计 + 置换区间。
这条限制应写进论文 Limitations。

**参考谱必须做 CPM 标准化**：4 类细胞的 RNA 含量差 4.9 倍，
不标准化会让体细胞比例高估 0.26，且偏倚方向恰好与待检验的组成位移同向。

---

## 内存与耗时

参考矩阵 37,241 × 34,633 稠密 float64 = **10.32 GB**，任何 `t()`/`scale()`/子集操作
都会再造一份 → 20.6 GB > R 的 `mem.maxVSize()` = 16 GiB，**必然失败**。

> 报错为 `vector memory limit of 16.0 Gb reached`，**看起来像"内存不够"，
> 很容易被误判为需要换机器**——实际只是方法错误。改用流式聚合后峰值仅 **47.7 MB**。

| 口径 | 量 |
|---|---|
| 机器计算合计 | 40–60 min（全局峰值 < 1.5 GB） |
| 人工核对与写作 | 8–15 h |
| 建议切分 | 7 个会话 × 2 h ≈ 3.5 周 |

最耗时的人工步骤是 **WT/KO 分组核对**：分组搞反不会报错，但结论完全反向。

---

## 已知限制

- **检验力**：见上。M = 6–8 无法分辨细粒度的解释比例。
- **参考系失配**：GSE112393 是成年小鼠，部分 bulk 是幼鼠（P0、P21）
  或不同品系 → 偏倚是系统性的，不是随机噪声。
- **校准的循环性**：用 GSE112393 合成伪 bulk 再用它去卷积必然成功，
  这不验证参考系匹配真实组织。必须补跨范式校准。
- **阻滞节点标注**：并非全部来自文献明确陈述，部分是按候选基因推断的。

---

## 引用

- Green CD, et al. *A Comprehensive Roadmap of Murine Spermatogenesis Defined by
  Single-Cell RNA-Seq.* Dev Cell. 2018. [doi:10.1016/j.devcel.2018.07.025](https://doi.org/10.1016/j.devcel.2018.07.025)
- Legendre P, Fortin MJ. *Comparison of the Mantel test and alternative approaches
  for detecting complex multivariate relationships in the spatial analysis of
  genetic data.* Mol Ecol Resour. 2010. [doi:10.1111/j.1755-0998.2010.02866.x](https://doi.org/10.1111/j.1755-0998.2010.02866.x)
- Lichstein JW. *Multiple regression on distance matrices: a multivariate spatial
  analysis tool.* Plant Ecol. 2007. [doi:10.1007/s11258-006-9126-3](https://doi.org/10.1007/s11258-006-9126-3)
- Newman AM, et al. *Robust enumeration of cell subsets from tissue expression
  profiles.* Nat Methods. 2015. [doi:10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337)
- Zhang J, et al. *Microarray meta-analysis identifies candidate genes for human
  spermatogenic arrest.* Andrologia. 2019;51(8):e13301. [doi:10.1111/and.13301](https://doi.org/10.1111/and.13301)（PMID 31037746）

完整引用见各文档文末。

---

## 许可

[MIT](LICENSE)
