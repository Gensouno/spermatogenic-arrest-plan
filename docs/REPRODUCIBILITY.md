# 复现、环境与验证范围

## 推荐路径

```bash
bash envs/setup_analysis.sh
bash envs/python.sh scripts/data_assets.py --download
bash envs/python.sh scripts/run_analysis.py --list
bash envs/python.sh scripts/run_analysis.py
```

需要 Python≥3.12、R≥4.4 和 curl；实测为 Python 3.14.7 / R 4.6.1、Apple Silicon。Python 顶层版本固定在 `envs/requirements-analysis.txt`。R 与 Bioconductor 必须匹配，安装脚本只补缺失包，不保证在不同 R 上自动得到同一注释版本。R 依赖包括 data.table、edgeR、limma、AnnotationDbi、org.Mm.eg.db、affy、mouse4302.db、mouse4302cdf、statmod。

依赖放在 `.venv/`、`.Rlib/`；不移动系统 Python/R，不使用 sudo。已有环境运行 `bash envs/setup_analysis.sh --check`。本机发布准备时这两个目录通过 Git 忽略的符号链接复用工作目录依赖，链接不提交到仓库；新用户应建立自己的环境，不依赖本机链接。

## 阶段与记录

19 个阶段依次为 counts、旧注释、扩展 RNA 导入与 DE、芯片导入与 DE、样本元数据、图谱聚合与校准、混合模型审计、图谱元数据、跨队列效应、组成敏感性、组成基因留出预测、稳健性与汇总、外部分选检查、细分阶段检查、代码测试。

`scripts/run_analysis.py` 保存 `logs/analysis/<UTC时间>/run.json`、每步输出及 Python 包列表；R 分析保存 sessionInfo。文件哈希追踪代码、分组表与输入清单。失败即停止，不把缺少队列静默忽略。`--from-stage` 仅用于已有前序产物的续跑，不自动重建被跳过的依赖；`--skip-asset-check` 只用于调试，并写入运行记录。

## 新目录复现与快照

```bash
# 先完成工作目录的全流程，再检查没有旧产物的情况
bash envs/python.sh scripts/check_clean_reproduction.py

# 完成后导出带校验值的小型结果，再检查
bash envs/python.sh scripts/export_snapshot.py
bash envs/python.sh scripts/export_snapshot.py --verify
```

新目录检查仅复制代码、分组表、冻结清单及 22 份输入，不复制旧结果或中间矩阵。重用已安装依赖，然后比较 18 张关键表：八张 DE 表、共同基因相关、敏感性表、组成预测与完整协方差、外部及细分参考诊断、混合模型审计两表。比较相对容差 1e-10、绝对容差 1e-12，并检查表结构和字符串。临时目录保留以便核查，路径在本地 `results/reproducibility/clean_run.json`；不将本机临时路径发布到结果快照。

快照 `manifest.json` 记录逐文件大小与 SHA-256、输入清单和分组表哈希、源代码哈希。`git_base_commit` 仅标识准备时的基线提交；实际准备中的未提交代码以 `source_sha256` 为准，避免把旧基线提交冒称为当前结果版本。更新分析代码后需重跑分析并重新导出快照，不可只修改图而保留旧清单。仅打包／启动脚本变动时，导出器会先确认原运行记录绑定的分析源码、输入说明和既有复现记录仍有效，再更新打包清单；这不代表重新运行了全部基线分析。

2026-10-23 新增本地可选代理读取后，已按上述打包流程更新清单，并逐文件确认原 48 项基线成果与前一发布版字节一致。可选数据扩展及其独立新目录检查另见 [数据扩展记录](DATA_EXPANSION.md)，不混用两次复现的范围。

快照不含完整逐基因矩阵。`validation_summary.json` 保存新目录比较的范围，并明确 `cell_number_fractions_validated=false` 和 `causal_composition_percentages_validated=false`。不能用“程序运行成功”替代生物学有效性。

2026-10-22 的发布准备中，正式目录的 19 阶段基线运行及无旧产物的新目录重跑均已完成，18 张关键结果表一致，9 项代码测试通过。详见快照的 `validation_summary.json`、环境 sessionInfo 和清单指纹。这是已安装依赖下的重复性证据，不是细胞比例通过外部验证的证据。

发布前可先对准备提交的文件执行 `bash envs/python.sh scripts/check_release.py`：核对快照、当前代码和输入说明指纹、Git 索引中的数据/环境排除、常见凭证模式及当前文档本地链接。历史副本中的相对路径保留原样，不作为当前运行链接检查目标。此扫描不能替代对所有潜在敏感信息的人工审阅。

## 尚未验证的范围

- 尚未在空依赖环境重新安装，也未验证其他系统和 R/Bioconductor 组合。
- 22 份输入已逐文件校验，但没有删除原数据后测试全部首次网络下载；curl 恢复下载也不是所有服务器都同等支持。
- 新目录运行复用现有包，所以不排除依赖版本变化导致差异。
- 代码测试数量有限；不等于全部实现已被形式化证明。没有原始表时，源 XLSX 测试会明确跳过，不算完成了数据验证。
- 旧大环境脚本、编号历史分析脚本与搜索发现脚本不属于当前默认路径。
