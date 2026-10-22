# 核查后结果快照

此目录只收录可直接查看的小型结果表及代表图。原始输入、完整逐基因差异表、完整中间矩阵和依赖环境不入库；可按根目录复现说明重算。快照不是完整实验数据的替代品。

入口：`cohort_summary.csv` 为八队列样本量和差异分析数量，`all_cohort_comparison/` 为八队列效应相关与 PCA，`robustness/` 为 63 个敏感性变体，其他目录为组成拟合校准及外部失配检查。

`manifest.json` 为逐文件 SHA-256、数据清单与分组表哈希、分析代码哈希；运行 `bash envs/python.sh scripts/export_snapshot.py --verify` 可检查文件完整性。`validation_summary.json` 区分新目录复现通过与科学假设成立。统计量和表格字段的解释见 [当前结果](../docs/RESULTS.md) 与 [方法](../docs/METHODS.md)。

拟合系数不能当作经过验证的细胞数比例，不报告因果组成贡献百分比。原始数据的作者、出处与许可不因本仓库的软件 MIT 许可而改变。除公开原表中明确标注的格式异常外，没有根据本次分析删掉任何样本。
