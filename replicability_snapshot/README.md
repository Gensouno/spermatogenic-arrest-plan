# 方向重复性与候选适用范围快照

2026-10-26 完成的探索性分支。五个计数队列、40 样本、四个遗传单位、9,675 共同受检基因；两个 Rfx2 年龄不重复计票。完整报告见 [REPLICABILITY.md](../docs/REPLICABILITY.md)。

- 正式结果：38,700 个方向／重复层级假设联合 BY 校正，发现 **0**；BH 敏感性同为 0。
- 描述性结果：105 个基因在至少三个单位有同向 1.5 倍变化，36 个下调基因在三个固定单位同时保留完整拟合及所有逐样本删除条件。这是后续优先级，没有筛选子表的 FDR 保证。
- 两来源阶段背景：36 个优先候选有 14 个大类一致、19 个不一致、3 个不能比较。不能把它们整体当作已证实的减数分裂机制。
- 六个既有额外对比的事后情境检查保留全部 105 候选和缺失，未增加主检验票数，也未重拟合外部研究。

| 文件 | 内容 |
|---|---|
| `directional_replicability_tests.csv.gz`、`replicability_summary.csv` | 两个完整宇宙的方向检验与校正 |
| `primary_BY_findings.csv` | 正式发现表；本次仅表头 |
| `primary_gene_summary.csv.gz`、`*_unit_evidence.csv.gz` | 全基因单位效应、P 值及稳定性 |
| `descriptive_1p5fold_candidates.csv.gz`、`candidate_stage_context.csv` | 完整 105 候选，不能当成显著表 |
| `candidate_unit_support.csv`、`stable_descriptive_priorities.csv` | 支持的固定单位、36 个优先候选及全部注释 |
| `all_primary_gene_stage_context.csv.gz`、`stage_*counts.csv`、`stage_context_summary.csv` | 全部主宇宙的阶段背景、来源分歧和分母 |
| `external_candidate_context.csv`、`external_context_*.csv` | 六个既有情境的逐基因覆盖、方向和留一范围 |
| `source_statistics_audit.csv`、`LOO_recomputation_audit.csv`、`R_adjustment_audit.csv` | 原拟合与40次计数留一对账、R/Python校正交叉核对 |
| `*design.json`、`runtime_*`、`reproduction.json` | 方法边界、实测版本、新目录复跑范围 |
| `manifest.json` | 文件、来源代码闭包和规格 SHA-256 |

新目录复跑比较 28 张表（含两张本地大表）、4 份 JSON 与3张图；83项测试通过。结果是计算复现与文件完整性证据，不是生物验证或空依赖安装验证。运行入口为 `scripts/run_replicability.py`；检查为 `scripts/check_replicability.py --inputs`。
