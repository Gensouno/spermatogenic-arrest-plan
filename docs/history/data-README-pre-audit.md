# 历史副本：不作为当前结果或运行入口

> 下文原样保留核查前的说明，其中部分结论已被推翻。文内相对路径按原文件位置解释。当前入口见仓库根 README.md 及 docs/HISTORY.md。

# 数据说明

本目录只放**清单与分组表**，原始表达数据**不入库**（体积大、且有 GEO 使用条款）。
按下面的步骤把数据取回来即可复现全流程。

## 本目录已有的文件

| 文件 | 内容 |
|---|---|
| `inventory.csv` | 220 行 GEO 自动盘点（检索 → 元数据 → 可用性判定）。**⚠️ 含 15 处误标**，见 `docs/02-datasets.md` |
| `dataset_selection.csv` | 57 行 × 19 列，逐条核验后的选择结果：core 8 / extended 4 / rejected 30 / mislabel 15 |
| `samples.tsv` | 86 行逐样本分组表（`tissue_class` 是去卷积的入口开关） |
| `README.md` | 本文件 |

## 下载后应当出现的路径

```
data/
├── raw/
│   ├── ref/
│   │   ├── GSE112393_DGE.txt.gz                 81,261,540 B   去卷积参考系（必须）
│   │   └── GSE112393_PerCellAttributes.txt.gz    1,495,253 B   细胞类型注释（必须）
│   ├── GSE236388/GSE236388_XAP5_counts.txt.gz      285,621 B   正向对照用（D16 睾丸臂）
│   ├── GSE283904/GSE283904_Gene_count.txt.gz       347,406 B
│   ├── GSE68283/GSE68283_Raw_counts.txt.gz         389,507 B
│   ├── GSE68284/GSE68284_Raw_counts.txt.gz         474,203 B
│   └── <其它 GSE>_RAW.tar                       见下表
└── ref/
    ├── signature_GSE112393.csv.gz                491,724 B   细胞类型特征谱（未标准化 counts）
    └── signature_GSE112393_cpm.csv.gz            384,037 B   ← 去卷积实际输入（CPM）
```

`data/raw/` 与 `data/ref/` 已在 `.gitignore` 中。

## 最小可用集合（跑通正向对照只需要两个文件）

1. **GSE112393 参考系**（82.8 MB）
   ```bash
   python3 scripts/fetch_geo.py \
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112393/suppl/GSE112393_MergedAdultMouseST25_DGE.txt.gz" \
     data/raw/ref/GSE112393_DGE.txt.gz

   python3 scripts/fetch_geo.py \
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112393/suppl/GSE112393_MergedAdultMouseST25_PerCellAttributes.txt.gz" \
     data/raw/ref/GSE112393_PerCellAttributes.txt.gz
   ```

2. **一个 bulk 数据集**（这里用 GSE236388 的 D16 睾丸臂）
   ```bash
   python3 scripts/fetch_geo.py \
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236388/suppl/GSE236388_XAP5_counts.txt.gz" \
     data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz
   ```

⚠️ **GSE236388 有两个 counts 文件，只有一个是睾丸**：

- `GSE236388_XAP5_counts.txt.gz` → **D16 睾丸**，KO 组标为 `XAP5 cKO`。**用这个。**
- `GSE236388_XAP5L_counts.txt.gz` → **D56 附睾尾精子**，不是睾丸，不能做组成分析。

3. 从参考系生成 CPM 特征谱（`data/ref/` 下的两个文件）
   ```bash
   MPLCONFIGDIR="$PWD/.mplcache" .venv/bin/python scripts/build_signature.py \
     --dge  data/raw/ref/GSE112393_DGE.txt.gz \
     --anno data/raw/ref/GSE112393_PerCellAttributes.txt.gz \
     --out  data/ref/signature_GSE112393.csv.gz \
     --group-out data/ref/cluster_collapse.csv

   MPLCONFIGDIR="$PWD/.mplcache" .venv/bin/python code/04_normalize_signature.py
   ```
   特征谱构建约 12 min；**峰值内存 ~48 MB**（流式，不会碰到 16 GB 上限）。

## 全部 8 个 CORE 数据集的体积

| GSE | 补充文件 | 体积 |
|---|---|---|
| GSE112393 | `*_DGE.txt.gz` + `*_PerCellAttributes.txt.gz` | 82.8 MB |
| GSE38096 | `GSE38096_RAW.tar` | 83.3 MB |
| GSE42601 | `GSE42601_RAW.tar` | 24.1 MB |
| GSE294209 | `GSE294209_RAW.tar` + 作者 DE 表 | 12.2 MB |
| GSE236388 | `GSE236388_XAP5_counts.txt.gz` | 0.3 MB |
| GSE101511 | `GSE101511_RAW.tar` | 1.9 MB |
| GSE283904 | `GSE283904_Gene_count.txt.gz` | 0.3 MB |
| GSE68283 | `GSE68283_Raw_counts.txt.gz` | 0.4 MB |
| GSE68284 | `GSE68284_Raw_counts.txt.gz` | 0.5 MB |

合计约 **205 MB**。

## 关于下载速度

NCBI GEO 的 FTP 对单连接限速很严（实测 **15–35 KB/s**）。`scripts/fetch_geo.py`
用 8 路并发 range 请求把聚合速率提到约 **126 KB/s**（约 6 倍），
并支持**块级断点续传**：中断后重跑只补缺失的块，全部就位后才拼接，
拼接前用 gzip CRC 校验，校验失败拒绝写入。所以按上面的命令反复重跑是安全的。

> 实测：81 MB 的 GSE112393 DGE 用 6–8 并发约 11–15 min。
> 若某轮遇到 NCBI 的 HTTP 503，脚本会退避重试，不需要人工干预。

## 数据来源与许可

所有数据来自 [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)，各数据集的使用请遵循
其原始发表论文与 GEO 的条款。本仓库不再分发这些文件。
