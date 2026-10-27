# 细线期参考边界与来源对账快照

研究判断、方法和局限见 [完整记录](../docs/REFERENCE_LIMITS.md)。本分支保留负结果，不改变原七组快照和八队列基线。

- `chen_*`：作者阶段与完整文库清单、逐细胞 QC、明确 ID 别名、来源阶段／文库分组；`chen_unresolved_cells.csv` 为空表，表示所有 1,136 个选中细胞已对账。
- `technical_*`、`gene_union_exon_lengths.csv.gz`：基因集合、权重和通用长度敏感性；保存所有类型，不能直接比较跨基因集合 RMSE。
- `training_marker_panel.csv.gz`：固定、仅由参考来源选出的阶段区分基因及排序证据，不依据外部目标标签。
- `transfer_*`、`prototype_evidence.csv`：五个候选的跨研究表达诊断，两个原型权重视角、两种基因权重及全基因／训练面板，均非独立重复。
- `source_library_*`、`source_holdout_controls.csv`：处理文库留出与外部敏感性，明确低细胞数组和精原细胞假阳性。
- `legacy_sorted_*`：原二十份分选资料的独立分支复查，preL 不参与大类评分。
- `decision_summary.json`：不把候选直接提升为已校准参考的证据摘要。
- `reproduction.json`、`manifest.json`：实际重跑边界和源代码／规格／结果指纹。

原始 UMI、补充材料 ZIP、完整阶段表达矩阵和候选参考矩阵在 Git 忽略目录中，不上传这些大文件。所有拟合权重均是表达兼容性诊断，不能作为已验证的体内细胞数比例。
