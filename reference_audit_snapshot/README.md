# 固定参考转移诊断快照（2026-10-24）

GSE121904 的 8 个发育期／富集捕获文库，生成全部导出／源代码 QC 与 RNA 求和／等细胞均值的 32 种目标谱；加上前轮 20 份分选 WT，合计 52 个诊断谱，不是 52 个动物或全睾丸重复。

`library_QC.csv`、`author_sample_table.csv` 和 `source_audit.json` 记录来源核对和标签缺口。作者时期分组不是逐细胞类型；没有提供新的独立年龄匹配标签真值。

`holdout_metrics.csv`、`fine_vs_broad.csv` 检查表达预测；`fold_weights.csv`、`weight_ranges.csv` 和 `design_conditions.csv` 检查系数稳定性。`sorted_class_matches.csv` 仅评分原已指定的 18 份谱，preL 不强制分类。不能把这些指标改写为细胞数比例准确率或独立实验显著性。

`reproduction.json` 记录新目录复跑范围；`manifest.json` 记录结果和源码指纹。详见 [方法、结果与限制](../docs/REFERENCE_TRANSFER.md)。上游文件保留原许可，作者 GPL 分析函数未复制进本项目。
