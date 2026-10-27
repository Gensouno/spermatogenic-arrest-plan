# 三个有作者标签文库的固定参考诊断

本目录只保存小型核查结果和图表；原始标注、矩阵字节窗口和完整中间表在 Git 忽略目录。

三个完整文库核对 8,791 个细胞，38 个至少 20 细胞的作者类型组形成 76 个表达处理谱。另有 EmptyDrops 元数据核查及 266 个作者偶线期重叠细胞的独立标注敏感性分析；未获取完整 EmptyDrops 表达，Leptotene 重叠仅 13 个、不评分。`reproduction.json` 记录 14 表、4 JSON 的新目录复跑范围。

范围预先限定为 E-MTAB-6946 的成年 B6、P5、P15 各一个完整正常文库；不是完整图谱。使用作者给出的细胞标签，不依据拟合结果重新命名细胞。RNA 求和与等细胞均值是同一批细胞的两种表达目标，不是独立重复。拟合权重不作为真实组织细胞数比例。

运行入口：`bash envs/python.sh scripts/run_labeled_reference.py`。离线检查：`bash envs/python.sh scripts/check_labeled_snapshot.py`，加 `--inputs` 同时核查缓存与冻结参考。新目录复跑：`bash envs/python.sh scripts/check_labeled_reproduction.py`。

科学解释和边界见 [说明](../docs/LABELED_REFERENCE.md)。`manifest.json` 验证内容和源代码／输入说明指纹，不等于科学有效性验证。
