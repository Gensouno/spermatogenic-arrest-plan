[English](README.md) | **简体中文**

# 生精阻滞跨模型转录组分析

比较小鼠全睾丸转录组中的共性与异质性，并检验单细胞参考能否支持细胞组成解释。项目使用公开数据再分析，不包含自行开展的动物或测序实验。

当前是经过数据和方法核查的探索性分析：8 个队列、56 个主分析样本、7 种遗传扰动。没有将“方向符合预期”当作比例准确性验证，也没有把残差直接称为共享分子机制。

## 主要结果

- 部分模型表达效应呈稳定正相关，也有明显反向关系；不支持直接把所有阻滞模型归为同一个转录组反应。
- 逐个删除 56 个样本，加上 TPM 伪计数及探针汇总变化，共检查 63 个变体。几个较强关联保持方向，Hnrnpgt 的弱关联不稳定；没有据此删除样本。
- 单细胞参考在独立分选资料的早期减数分裂群体上存在失配；细分至 12 个生殖状态后仍未充分解决。因而不报告已验证的细胞数比例或因果组成贡献百分比。

![八队列共同受检基因的 KO−WT 效应相关](results_snapshot/all_cohort_comparison/eight_cohort_effects.png)

图中为 8,099 个共同受检基因的 Spearman 相关，不是对独立生物重复的显著性检验。TPM 与芯片作为敏感性分支，不能和计数效应完全等量解释。完整说明见 [结果与边界](docs/RESULTS.md)，精简表格见 [结果快照](results_snapshot/README.md)。

## 可选数据扩展（2026-10-23）

原八队列基线保持不变。可选探索分支定向核查了 32 个 GEO 系列，接入 Dcaf2 P9、Pramef12 P2/P7（杂合子对照）和 Msl3 P45（表型阴性比较），新增 24 份全睾丸样本，另检查 20 份分选 WT 表达谱。不同角色不混为“新增阻滞病例”。方法、文件身份问题及复现范围见 [数据扩展记录](docs/DATA_EXPANSION.md)，结果在 `expansion_snapshot/`。

## 数据与分析范围（冻结的八队列基线）

2026-10-26 另完成 [Setx 发育序列与 Eef1g 分选细胞的新来源核查](docs/NEW_MODEL_CONTEXTS.md)：两条GEO查询返回61个去重系列，核查其中12个来源，下载并单独分析40份Setx全睾丸TPM及6份Eef1g分选L/Z计数谱。固定36个下调优先候选在Setx P14/P16/P18分别24、19、33个经留一仍保持同向1.5倍，但不能排除细胞构成；没有增加主检验票数。新目录复跑14张表、2份JSON和2张图一致；结果为 `new_models_snapshot/`，入口为 `scripts/run_new_models.py`。

2026-10-26 完成 [跨遗传模型方向重复性与候选适用范围](docs/REPLICABILITY.md)：五个计数队列的 40 份样本按四个遗传单位检验，Rfx2 两年龄不重复计票。9,675 基因的完整方向假设族经 BY／BH 校正均为零发现；单列 105 个描述性候选及 36 个逐样本删除后稳定的下调优先候选，没有称作已证实共享机制。两个单细胞来源的阶段分歧和六个既有额外对比的事后覆盖检查均完整保留。新目录复跑 28 张表、4 份 JSON、3 张图一致；结果为 `replicability_snapshot/`，入口为 `scripts/run_replicability.py`。

2026-10-25 完成 [细线期参考边界与独立 UMI 对账](docs/REFERENCE_LIMITS.md)：从 GSE107644 原始矩阵及作者 S1/S2 核对 1,136 个阶段明确细胞，补充基因集合、权重、长度敏感性和来源文库留出。简单增加 L/Z 参考不能稳定纠正外部失配，且会把部分精原细胞主要拟合为减数分裂；回查原二十份分选谱也未解决 Cnot4 LZ。结果单列 `reference_limits_snapshot/`，没有把候选提升为已校准细胞比例模型。检查入口为 `scripts/check_reference_limits.py`。

2026-10-25 继续完成 [独立细线期来源诊断](docs/LEPTOTENE_FOLLOWUP.md)：核查四个候选研究，接入 GSE155967 阶段分开的六份 RNA 计数和 GSE83264 的四份混合 LZ／粗线期对照。两份明确细线期在细分参考下仍主要拟合为精原细胞，两份偶线期则转为减数分裂大类；保留同步化、年龄、分选、基因型记录及文库技术限制。PreL 不参与大类标签评分，没有据此训练新参考或改动八队列基线。结果为 `leptotene_snapshot/`，检查入口为 `scripts/check_leptotene.py`。

2026-10-25 补齐 [完整 P15 EmptyDrops 早期细胞诊断](docs/EMPTYDROPS_P15.md)：逐细胞核对 13,168 个作者标注细胞。细分参考对全部 938 个偶线期、其中 672 个新增细胞的改善延续，但 907 个细线期仍主要拟合为精原细胞，不能宣布早期减数分裂失配整体解决。同条形码、同标签和同基因对照中，两版计数的汇总与拟合结果一致；新目录复跑 14 张表与两份设计／来源 JSON 一致。新结果单列 `emptydrops_p15_snapshot/`，不改变八队列基线；核查入口为 `scripts/check_emptydrops_p15.py`。

2026-10-24 补充 [作者标签明确的三文库参考适配核查](docs/LABELED_REFERENCE.md)：完整核对成年 B6、P5、P15 共 8,791 个细胞，诊断作者标注的 38 个类型组；并定位到明确的细线／偶线期作者标签。已对有缓存表达的 266 个偶线期重叠细胞作独立敏感性分析，但大多数低复杂度细胞表达尚未获取，不能宣布 LZ 失配或细胞比例准确性已解决。快照为 `labeled_reference_snapshot/`，检查入口为 `scripts/check_labeled_snapshot.py`。

2026-10-24 另完成 [发育期来源核查与固定参考转移诊断](docs/REFERENCE_TRANSFER.md)：检查 8 个捕获文库及 15,882 个细胞，对 20 份已有分选谱与 32 种文库处理组合做留出基因预测。细分参考的预测误差有所改善，但未据此宣布细胞比例准确；仍缺逐细胞作者标签。统一离线完整性检查入口为 `scripts/check_project.py`。

2026-10-24 另完成 [Lrh1 早期对比及 Ppp2ca 文件核查](docs/FOLLOWUP_LRH1.md)：16 份 Lrh1 样本分为两组独立比较，处理作者批次并保留方法敏感性；Ppp2ca 的转录本 FPKM 仅核查、不混入计数主分析。这部分为单独的后续探索，快照在 `followup_snapshot/`，不改变以下八队列范围。

| 队列 | 遗传扰动 | KO / WT | 数据类型 |
|---|---|---:|---|
| GSE236388 | XAP5，D16 睾丸臂 | 3 / 3 | 基因 counts |
| GSE283904 | N6amt1 | 4 / 4 | 基因 counts |
| GSE68283 | Rfx2，P29–30 | 4 / 5 | 基因 counts |
| GSE68284 | Rfx2，P21 | 6 / 5 | 基因 counts |
| GSE101511 | Hnrnpgt，18 dpp | 3 / 3 | 旧转录本特征 counts |
| GSE294209 | Rnu11，P14 | 3 / 3 | 作者提供的 TPM |
| GSE38096 | Sin3a，P0 | 2 / 2 | Agilent 单 Cy3 通道 |
| GSE42601 | Dazap1，P21 | 3 / 3 | Affymetrix mouse4302 |

Rfx2 两年龄队列来自同一研究，不视作两种独立遗传模型。Sin3a 两个 HET 不纳入主对比；XAP5L 的 D56 附睾尾精子不纳入睾丸分析。N6amt1 年龄未由当前官方样本记录确认，不沿用旧表中未经确认的 adult 标签作年龄分析。阻滞分期并非每队列都已获得原始组织学证据，不据旧节点标签作分期因果推断。

各队列独立处理，合并的是基因效应而不是原始表达矩阵。组成探索只用于五个 counts 队列。22 份输入的官方 URL、大小与 SHA-256 见 [资产清单](data/analysis_assets.json)；来源和分组说明见 [数据说明](data/README.md)。

## 复现

需要已安装的 Python、R 和 curl。当前实测 Python 3.14.7 / R 4.6.1；不同 R/Bioconductor 注释版本可能改变基因映射。依赖安装不需要移动系统解释器或使用 sudo。

```bash
git clone https://github.com/Gensouno/spermatogenic-arrest-plan.git
cd spermatogenic-arrest-plan

# 首次安装；已有环境可加 --check 仅检查
bash envs/setup_analysis.sh

# 下载缺失输入并检查冻结校验值
bash envs/python.sh scripts/data_assets.py --download

# 19 个阶段，从输入重算
bash envs/python.sh scripts/run_analysis.py

# 检查随仓库提供的结果快照是否完整
bash envs/python.sh scripts/export_snapshot.py --verify
```

已有输入时用 `scripts/data_assets.py --verify` 检查；也可通过 `--import-from /path/to/existing-project` 从另一目录复制经校验的输入。已有输入若校验不符会停止，不自动覆盖。

### 可选下载代理

需要代理时，将 `envs/download-proxy.example` 复制为项目根目录的 `.download-proxy.local`，把其中 URL 改为本机 HTTP(S) 代理地址。该文件已被 Git 忽略，不提交代理地址或凭据。通过 `envs/python.sh`／`envs/rscript.sh` 运行的下载及依赖安装会继承它；代理不可用时应修复代理，不自动改为直连。未配置此文件时保留调用者的代理环境。直接运行系统 Python/R 不会读取此项目配置。

本机另做了没有旧结果和中间文件的新目录重跑，关键结果表一致；此检查复用已安装的依赖，**不是**另一台机器或空依赖环境安装验证，也没有在删除原数据后测试全部首次网络下载。耗时和详细边界见 [复现记录](docs/REPRODUCIBILITY.md)。

## 项目导航

- [Setx 年龄交互分析](docs/SETX_INTERACTION.md)：联合年龄／基因型模型、来源编号和条码限制、敏感性与复现。
- [Setx 来源与长度核查](docs/SETX_SOURCE_AND_LENGTH.md)：数据库身份、等位基因背景和长度校正计数敏感性（英文报告）。
- [Setx 逐样本稳定性](docs/SETX_STABILITY.md)：计数模型留一检查、同基因集合诊断及模型依赖（英文报告）。
- [新增全睾丸背景](docs/INDEPENDENT_CONTEXTS.md)：Rif1、Hormad1/2、Mael 和 Morc1 来源核查与独立计数分析；原有主要检验集合保持不变（英文报告）。
- [英文论文同步工作流](docs/PAPER_WORKFLOW.md)：持续写作、证据登记、检查、Word导出和私稿发布边界。
- [METHODS](docs/METHODS.md)：分平台处理、基因集合、敏感性与组成检查。
- [RESULTS](docs/RESULTS.md)：主要数字、负结果及解释边界。
- [REPRODUCIBILITY](docs/REPRODUCIBILITY.md)：运行顺序、环境与验证范围。
- [面试讲解提纲](docs/INTERVIEW.md)：帮助理解自己使用的数据、代码和局限。
- [HISTORY](docs/HISTORY.md)：早期 AI 草稿、被推翻的假设与保留原因。

`code/`、`scripts/` 为当前实现；`results_snapshot/` 为带校验值的精简成果；`results/` 为重算生成的完整产物，不入库。原始数据、依赖环境和大矩阵不提交。

`PLAN.md`、编号 `docs/01`–`04`、`notes/` 及三个编号旧分析脚本是历史材料，不是当前执行方案。

## 来源与许可

单细胞参考为 [Green et al., Dev Cell 2018](https://pubmed.ncbi.nlm.nih.gov/30146481/)，外部分选资料为 [da Cruz et al., BMC Genomics 2016](https://link.springer.com/article/10.1186/s12864-016-2618-1)。各队列出处见资产清单与逐样本 GEO 证据。

原创代码使用 [MIT](LICENSE) 许可；上游数据、原论文与补充材料保留其原有许可和引用要求，不因本仓库软件许可而改变。
