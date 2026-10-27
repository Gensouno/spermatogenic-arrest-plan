# 独立细线期来源诊断快照

本轮单列两个新研究的 10 份 RNA 表达谱：GSE155967 的 PreL/L/Z 各两份，GSE83264 的混合 LZ/粗线期各两份。PreL 不参与大类标签评分，其余八份仅作固定参考适配诊断，不计为新增遗传阻滞队列。方法与局限见 [研究记录](../docs/LEPTOTENE_FOLLOWUP.md)。

- `candidate_review.csv`：四个候选来源的纳入、延后及原因。
- `sample_evidence.csv`：逐样本来源、阶段、制备、基因组、计数单位及表达量覆盖。
- `mapping_audit.json`：Ensembl 映射损失和排除的 HTSeq 统计行。
- `replicate_concordance.csv`：每来源／阶段两份样本的描述性表达相关。
- `stage_comparison.csv`、`full_fit_weights.csv`、`fold_weights.csv`、`holdout_metrics.csv`：全部处理结果，包括负结果和未评分 PreL。
- `gene_design.csv.gz`、`design_conditions.csv`、`design.json`：固定基因、五折、权重、条件数和原文来源证据。
- `stage_diagnostics.png`：RNA 权重与表达误差对照。
- `reproduction.json`、`manifest.json`：实际重跑范围及文件／源代码／输入规格指纹。

原始计数、来源 SOFT、论文 XML 与完整共有基因计数矩阵位于忽略目录，不上传大文件。两来源技术和年龄不同，不能用跨来源误差高低推断生物学质量；RNA 拟合权重不能解释成细胞数比例。
