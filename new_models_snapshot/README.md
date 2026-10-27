# 新来源与候选表达背景

2026-10-26：定向查询GEO得到61个去重系列，完成12个来源的样本级核查；下载并分析Setx五年龄40份全睾丸TPM谱、Eef1g P12六份分选L/Z计数谱。报告见 [NEW_MODEL_CONTEXTS.md](../docs/NEW_MODEL_CONTEXTS.md)。

- Setx的40份基因RSEM表与40份piRNA表分开，保存估计计数、长度和TPM；没有将估计值四舍五入作为原始计数。
- 固定36个下调优先候选在Setx P14/P16/P18分别24、19、33个经留一后仍保持同向1.5倍；这是描述性情境延续，细胞构成仍可参与。
- Eef1g仅5/36个候选受检且唯一映射，不能据少量覆盖下结论。
- 两来源未加入四遗传单位主检验，没有改变正式零发现结论；动物独立性、合池、阶段比例和隐藏批次尚未充分核实。

`search_inventory.csv`、`query_audit.json`、`metadata_review.csv`、`new_assets.json`保留检索与下载身份；`sample_evidence.csv`和`Setx_archive_audit.csv`保留样本及成员对账；`candidate_context.csv`、`context_summary.csv`、`gene_coverage.csv`及图保留完整候选和缺失；`source_mapping_mass.csv`保留注释损失。大矩阵与完整条件模型统计在本地Git忽略目录。

新目录复跑14张表、2份JSON和2张图一致，87项全套测试通过。运行 `scripts/run_new_models.py`；检查 `scripts/check_new_models.py --inputs`。完整性／计算复现不是动物独立性、细胞比例或机制的生物验证。
