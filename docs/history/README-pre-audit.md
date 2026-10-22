# 历史副本：不作为当前结果或运行入口

> 下文原样保留核查前的说明，其中部分结论已被推翻。文内相对路径按原文件位置解释。当前入口见仓库根 README.md 及 docs/HISTORY.md。

# 生精阻滞模型的转录组共性分析

本项目分析小鼠生精阻滞（spermatogenic arrest）模型的全睾丸 bulk 转录组，关注一个容易被忽略的问题：不同模型中反复出现的“共同差异基因”，在多大程度上反映共享的分子调控，在多大程度上仅由细胞组成变化造成。

全睾丸 bulk RNA-seq 是混合组织信号。晚期生精细胞减少会使其高表达基因整体降低，同时提高体细胞在总体信号中的相对贡献。因此，跨模型的相似表达变化不能直接等同于共同机制。本项目通过细胞组成零模型对这一解释进行定量检验。

## 项目范围与当前状态

项目提供分析方案、方法学规格和可复现代码。当前已完成并报告了一个去卷积正向对照；跨模型共性分析及其最终生物学解释仍需以完整运行结果、稳健性分析和独立验证为依据。

| 分析阶段 | 目标 | 当前说明 |
|---|---|---|
| 跨模型差异表达比较 | 描述不同模型的 logFC 是否收敛 | 分析框架已定义 |
| 细胞组成零模型 | 估计仅由组成变化可产生的表达差异 | 方法规格已定义 |
| 残差信号分析 | 在模型假设下检验未被组成预测解释的信号 | 依赖前述分析的完整结果 |
| 去卷积正向对照 | 检查端到端流程是否产生方向合理的估计 | 已运行，见下文 |

上述状态不应被解读为对最终科学结论的验证。

## 已运行的去卷积正向对照

以真实数据进行了端到端去卷积测试：

```text
数据集：GSE236388（XAP5 cKO 睾丸臂，D16；3 KO / 3 WT）
参考系：GSE112393（成年小鼠睾丸；37,241 genes × 34,633 cells，聚合为 4 类细胞）
```

| 指标 | 结果 |
|---|---|
| 重建 R² | 0.680–0.742（预设阈值：0.6） |
| germ-cell 比例 ΔP | −0.194（KO 低于 WT） |
| somatic 比例 ΔP | +0.194（KO 高于 WT） |

结果与该模型中生殖细胞减少、体细胞相对比例升高的预期方向一致。这是对流程可运行性和方向一致性的正向证据，而非对细胞比例绝对准确性或参考谱适配性的充分验证。复现命令见“快速开始”。

## 快速开始

```bash
git clone <repo-url>
cd spermatogenic-arrest-plan

# 1) 环境自检（幂等，约 15 s）
bash envs/setup.sh

# 2) 复现已运行的正向对照（约 1 min）
#    需先下载两个数据文件，见 data/README.md
MPLCONFIGDIR="$PWD/.mplcache" .venv/bin/python code/06_validate_deconvolution.py \
    --sig  data/ref/signature_GSE112393_cpm.csv.gz \
    --bulk data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz \
    --wt-pattern WT --ko-pattern cKO
```

环境要求：macOS aarch64（Apple Silicon）、R ≥ 4.4、Python ≥ 3.10 和 16 GB 内存。`envs/setup.sh` 将 R 包安装至项目内 `.Rlib`，Python 包安装至 `.venv`，不写入系统目录；如全局 R 库可写，也可使用系统安装。

GEO 下载受 NCBI 限速影响（单流约 15–35 KB/s）。仓库提供支持断点续传的 8 并发分块下载器 `scripts/fetch_geo.py`，实测约可加速 6 倍。原始数据不纳入仓库，下载说明见 [`data/README.md`](data/README.md)。

## 分析设计

分析分为三个相互衔接的阶段：

1. 比较不同模型的差异表达方向，描述跨模型 logFC 的收敛程度。
2. 构建仅包含细胞组成变化的零模型，量化其对跨模型相似性的可解释部分。
3. 将观察到的变化与组成预测进行比较，研究在模型假设下仍未被解释的残差信号。

零模型将细胞比例位移通过参考表达谱推进至 bulk 表达空间：

$$\hat c_m(g)=\log\!\big((S\,\bar p^{KO}_m)(g)+1\big)-\log\!\big((S\,\bar p^{WT}_m)(g)+1\big)$$

其中，$\hat c_m(g)$ 表示：若基因的细胞内表达不变、仅细胞组成改变，模型 $m$ 中基因 $g$ 预期呈现的 logFC。该模型估计的是给定参考谱、细胞类型划分及模型假设下由组成变化可解释的表达份额；它本身不能确定因果机制，也不能将残差信号直接解释为直接调控。

### 方法学要点

- 直接将 logFC 投影到 $S$ 的列空间并不能正确估计组成效应：即使使用完美参考谱，模拟中也仅得到 0.585（真值 1.0），系统性低估约 40%。原因是 log 为凹函数，组成效应不落在 $S$ 的列空间中。
- 朴素 Mantel 检验在“纯伪影”与“存在真实信号”的模拟情形下均给出 *p* = 0.158，区分能力不足；因此采用以模型对为单位的 MRM 回归。
- 置换单位应为“模型”，并采用模型错配。基因打乱法的零分布 sd = 0.007，过窄且易产生大量假阳性；模型错配法的 sd = 0.129。
- 当 M = 6–8 时，分析可用于判断组成变化是否为主要解释之一，但不足以可靠地区分“37%”与“44%”等细粒度解释比例。因此 `compShare` 应报告点估计和置换区间。
- 参考谱须进行 CPM 标准化。四类细胞的 RNA 含量相差 4.9 倍；不标准化时体细胞比例会被高估 0.26，且偏倚方向与待检验的组成位移一致。

数学规格、校准准则和不确定性清单见 [`docs/03-null-model.md`](docs/03-null-model.md)。

## 数据

分析纳入 8 个全睾丸生精阻滞模型；每个模型均按 GEO 原始信息逐样本核验。

| GSE | 模型 | 阻滞节点 | n（KO/WT） | 年龄 | 平台 |
|---|---|---|---|---|---|
| GSE38096 | Sin3a 生殖细胞特异 cKO | 精原细胞/起始 | 6（2/2 + 2 HET） | P0 | GPL13912 |
| GSE42601 | Dazap1 突变 | 精原细胞/起始 | 6（3/3） | P21 | GPL1261 |
| GSE294209 | Rnu11 cKO | 减数分裂前期 | 6（3/3） | P14 | GPL24247 |
| GSE236388 | XAP5 cKO（仅 D16 睾丸臂） | 减数分裂前期 | 6（3/3） | D16 | GPL24247 |
| GSE101511 | Hnrnpgt cKO | 减数分裂前期 | 6（3/3） | 18 dpp | GPL19057 |
| GSE283904 | N6amt1-SKO | 减数分裂前期 | 8（4/4） | 成年 | GPL24247 |
| GSE68283 | Rfx2 −/− | 减数分裂后/精子形成 | 9（4/5） | P29–30 | GPL17021 |
| GSE68284 | Rfx2 −/−（同论文 P21） | 减数分裂后/精子形成 | 11（6/5） | P21 | GPL17021 |

去卷积参考系为 **GSE112393**：成年小鼠睾丸单细胞图谱（Green et al., *Dev Cell*, 2018）。

> 去卷积仅适用于全睾丸样本。分选细胞群体的组成由分选步骤决定，并非待估计的生物学比例；`code/05_deconvolution.R` 对此设置了硬性 `stop()` 守卫。

公开数据集标签存在误标。`data/inventory.csv` 为自动盘点结果，其中有 **15 条**经逐条核验后确认存在问题，包括将睾丸类器官列为全睾丸、将附睾尾精子混入睾丸臂，以及样本数记录错误。分析应以 `data/dataset_selection.csv` 为准；纠错详情见 [`docs/02-datasets.md`](docs/02-datasets.md)。

## 目录结构

```text
.
├── PLAN.md                     总方案（建议从此开始）
├── docs/
│   ├── 01-environment.md       本机环境实测、内存峰值和已知限制
│   ├── 02-datasets.md          CORE-8 数据集筛选证据链与公开数据纠错
│   ├── 03-null-model.md        零模型数学规格、校准判据和不确定性清单
│   └── 04-pipeline-runbook.md  端到端执行规格（Step 0–9，含失败处理）
├── code/
│   ├── 04_normalize_signature.py     参考谱 CPM 标准化
│   ├── 05_deconvolution.R            NNLS 去卷积、硬门槛和 held-out R²
│   └── 06_validate_deconvolution.py  端到端验证器与正向对照
├── scripts/
│   ├── 01_geo_search.py … 03_curate.py  数据集盘点
│   ├── build_signature.py                流式构建细胞类型特征谱
│   └── fetch_geo.py                      并发分块下载器
├── envs/setup.sh               幂等环境搭建
├── data/                       数据集清单与逐样本分组表（不含原始表达数据）
└── notes/                      两版早期构想与数据集盘点全文
```

## 资源要求与执行注意事项

参考矩阵大小为 37,241 × 34,633。以稠密 float64 存储时约需 **10.32 GB**；`t()`、`scale()` 或子集操作可能产生额外副本，使内存需求升至约 20.6 GB，超过 R 的 `mem.maxVSize()` = 16 GiB，因而会失败。

典型报错为 `vector memory limit of 16.0 Gb reached`。这反映的是实现方式的内存复制问题，而非必然需要更换硬件。采用流式聚合后，峰值内存为 **47.7 MB**。

| 项目 | 估计 |
|---|---|
| 机器计算合计 | 40–60 min（全局峰值 < 1.5 GB） |
| 人工核对与写作 | 8–15 h |
| 建议切分 | 7 个会话 × 2 h，约 3.5 周 |

WT/KO 分组核对是关键人工步骤：分组颠倒通常不会触发程序错误，但会使结论方向完全相反。

## 局限性

- **统计检验力：** M = 6–8 的模型数限制了解释比例的精细推断能力。
- **参考系失配：** GSE112393 来自成年小鼠；部分 bulk 样本为幼鼠（P0、P21）或不同品系。该偏倚可能是系统性的，而非随机噪声。
- **校准循环性：** 以 GSE112393 合成伪 bulk 并再以其去卷积会自然获得较好结果，不能验证参考谱与真实组织的匹配；仍需跨范式校准。
- **阻滞节点标注：** 并非全部直接来自文献明确陈述；部分依据候选基因推断。

## 引用

- Green CD, et al. *A Comprehensive Roadmap of Murine Spermatogenesis Defined by Single-Cell RNA-Seq.* Dev Cell. 2018. [doi:10.1016/j.devcel.2018.07.025](https://doi.org/10.1016/j.devcel.2018.07.025)
- Legendre P, Fortin MJ. *Comparison of the Mantel test and alternative approaches for detecting complex multivariate relationships in the spatial analysis of genetic data.* Mol Ecol Resour. 2010. [doi:10.1111/j.1755-0998.2010.02866.x](https://doi.org/10.1111/j.1755-0998.2010.02866.x)
- Lichstein JW. *Multiple regression on distance matrices: a multivariate spatial analysis tool.* Plant Ecol. 2007. [doi:10.1007/s11258-006-9126-3](https://doi.org/10.1007/s11258-006-9126-3)
- Newman AM, et al. *Robust enumeration of cell subsets from tissue expression profiles.* Nat Methods. 2015. [doi:10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337)
- Zhang J, et al. *Microarray meta-analysis identifies candidate genes for human spermatogenic arrest.* Andrologia. 2019;51(8):e13301. [doi:10.1111/and.13301](https://doi.org/10.1111/and.13301)（PMID 31037746）

完整引用见各文档文末。

## 许可

[MIT](LICENSE)
