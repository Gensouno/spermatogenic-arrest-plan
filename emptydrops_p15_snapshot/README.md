# 完整 P15 EmptyDrops 诊断快照

来源、分析范围和解释边界见 [研究记录](../docs/EMPTYDROPS_P15.md)。本目录与原五个快照分开保存，不改变八队列基线。

- `group_inventory.csv`：包括作者 Outliers 与小于 20 细胞的排除组。
- `profile_evidence.csv`、`cell_complexity_summary.csv`：分支、细胞集合、作者类型、表达加权与复杂度证据。
- `early_meiotic_comparison.csv`：L/Z 及粗线期固定参考结果；不足阈值的重叠细线期不评分。
- `matched_count_branch_comparison.csv`：完全相同条形码、标签及基因宇宙的计数处理对照。
- `holdout_metrics.csv`、`full_fit_weights.csv`、`fold_weights.csv`、`stratum_summary.csv`、`design_conditions.csv`：全部合格组的诊断，不只展示焦点细胞。
- `gene_design.csv.gz`：所有本轮拟合共用的基因、固定折数和权重。
- `preparation_audit.json`、`design.json`：逐细胞核对与固定设计记录。
- `reproduction.json`：新目录实际重跑范围。
- `manifest.json`：结果和来源代码／规格指纹。

原矩阵窗口、逐细胞 QC 和类型表达矩阵在忽略目录中，不把原始大文件上传仓库。一个文库的细胞集合及两种加权视角不作为独立动物重复；RNA 拟合系数不作为已验证细胞数比例。
