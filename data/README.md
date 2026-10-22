# 数据与分组说明

默认分析只使用 `analysis_assets.json` 中冻结的 22 份输入，约 250.6 MiB。文件包括八个 bulk 队列、GSE112393 图谱及其元数据、Ensembl 79 注释、外部分选补充表。精子臂、未使用的备用队列和搜索发现文件不在默认输入清单中。

每项资产都有准确 URL、字节数及 SHA-256。GEO 补充文件 URL 与官方 SOFT 的 supplementary-file 字段核对；哈希代表本次获取的快照，不是 GEO 官方对内容正确性的认证。

```bash
bash envs/python.sh scripts/data_assets.py --download
bash envs/python.sh scripts/data_assets.py --verify
```

下载仅补齐缺失文件，校验通过才提升为正式输入；已有输入不符时停止，保留原文件。代理如有需要由调用者环境配置，仓库不写死本机代理地址。首次完整网络下载路径尚未在空目录验证。

`samples.tsv` 保留早期盘点的完整分组表，包含默认分析以外的候选样本；它不代表全部样本都进入分析。当前 56 个主分析样本以 `results_snapshot/import_qc/official_sample_evidence.csv` 为准，其中保存 assigned_group、原始 GSM 标题与来源信息。自动步骤核对样本 ID、标题及 Hnrnpgt 基因型；其余分组还依赖核查时对官方特征的人工解释，不把单纯标题匹配称为全部分组的自动证明。

重要边界：

- GSE236388 只选 XAP5 D16 睾丸，不选 XAP5L D56 附睾尾精子。
- GSE38096 主比较为 2 KO / 2 WT，两个 HET 保留但不纳入。
- GSE68283 与 GSE68284 是同一 Rfx2 研究的两个年龄队列。
- GSE294209 仅有每样本 TPM，不作为原始 counts。
- GSE101511 的 ENSMUST 旧特征使用固定 Ensembl 79 兼容映射；原始 GTF 发布版本未明确，不能称为完全复现作者注释。
- 旧表中 N6amt1 的 adult 年龄及部分“阻滞节点”缺少当前充分证据，不用于年龄或组织学分期因果分析。

`inventory.csv`、`dataset_selection.csv` 是历史搜索盘点，不是新的纳入证据；旧说明备份见 [历史数据说明](../docs/history/data-README-pre-audit.md)。

`data/raw/`、`data/ref/`、`data/bulk/`、`data/interim/` 仅在本地生成并被 Git 忽略。上游数据、原论文及补充表不适用本仓库原创代码的 MIT 许可；使用时遵守各来源要求并引用原作者。
