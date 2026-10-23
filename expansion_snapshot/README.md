# 可选扩展快照（2026-10-23）

此目录与 `results_snapshot/` 的八队列基线分开。新增分析为 DCAF2 P9、PRAMEF12 P2/P7（杂合子对照）和 Msl3 P45（表型阴性比较），共 24 份全睾丸样本；另检查三个研究的 20 份分选 WT 表达谱，不将这些谱当成额外全睾丸样本或已知比例混合物。

`cohort_summary.csv` 为样本角色与基因映射信息；`fixed_core_gene_correlations.csv` 将新增效应限制在原基线共同基因集合内，每行另报告实际共有基因数。`cross_cohort_ranges.csv` 只删除新增研究中的一个样本，原研究保持基线，范围不是置信区间。

`GSE232415_file_identity_audit.csv` 记录补充文件的内部标签不一致与重复向量；不以表达模式猜测修复标签。该系列全睾丸部分暂不纳入，单独检查过的四份分选 WT 谱只作参考兼容性诊断。未发现问题不等于证明所有实验身份准确。

完整筛选理由见 `data/expansion_review.json` 和 `docs/DATA_EXPANSION.md`。下载来源与 SHA-256 见 `data/expansion_assets.json`；本目录 `manifest.json` 保存结果和分析代码指纹。数据仍有对照、年龄、注释、技术平台和重复独立性的限制，不据此报告因果机制或可靠细胞比例。

`reproduction.json` 记录新目录复跑的 12 张表一致性检查。复跑使用缓存原始输入、冻结的基线效应／参考产物，并共享已安装依赖；不是空环境安装或全基线复跑。
