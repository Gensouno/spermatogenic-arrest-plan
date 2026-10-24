# 后续核查快照（2026-10-24）

两组 Lrh1 早期全睾丸对比，共 16 份样本；另核查 Ppp2ca 四个转录本 FPKM 文件，不纳入计数分析。原八队列与前轮四对比扩展不覆盖。

`cohort_summary.csv`、`Lrh1_samples.csv`、两个设计矩阵记录主比较及控制等位基因／作者批次；`batch_sensitivity.csv` 和 `method_sensitivity.csv` 区分预先保留的批次敏感性与显著性不一致后增加的方法诊断。`paper_direction_checks.csv` 不是整体准确率验证，不把方向一致写为四基因显著性复现。

`baseline_correlations.csv` 限定于旧 8,099 基因宇宙再取唯一映射受检交集；`baseline_ranges.csv` 只删除当前 Lrh1 样本，范围不是置信区间。`Ppp2ca_unit_audit.csv` 仅记录文件身份和单位。

`reproduction.json` 记录共享依赖、缓存原始输入及冻结基线效应的新目录表格比较，不是空环境安装或旧分析重跑。`manifest.json` 记录小型成果及源代码指纹。完整边界和命令见 [后续分析记录](../docs/FOLLOWUP_LRH1.md)，输入 URL 和校验值见 [资产清单](../data/followup_assets.json)。
