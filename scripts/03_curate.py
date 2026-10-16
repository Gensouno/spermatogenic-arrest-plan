#!/usr/bin/env python3
# 第三步：把 main / atlas 候选整理成可用的数据集清单
# 输入 data/geo_meta.tsv + data/geo_candidates.tsv
# 输出 data/inventory.csv 和 notes/数据集清单.md
#
#   python3 scripts/03_curate.py

import csv, os, re, sys
from collections import Counter, defaultdict

DATA = "data"
ROOT = "."

# ---------- 基因 -> 阻滞节点（依据文献常识，作为初判，需人工复核） ----------
GENE_NODE = {
    # 精原细胞/分化起始
    "Dmrt1": "精原细胞/起始", "Dazl": "精原细胞/起始", "Sohlh1": "精原细胞/起始",
    "Sohlh2": "精原细胞/起始", "Zbtb16": "精原细胞/起始", "Plzf": "精原细胞/起始",
    "Gfra1": "精原细胞/起始", "Nanos2": "精原细胞/起始", "Nanos3": "精原细胞/起始",
    "Neurog3": "精原细胞/起始", "Sall4": "精原细胞/起始", "Tfap2c": "精原细胞/起始",
    "Cxxc1": "精原细胞/起始", "Ezh1": "精原细胞/起始", "Chd8": "精原细胞/起始",
    "Dcaf2": "精原细胞/起始", "Mgat2": "精原细胞/起始", "Pramef12": "精原细胞/起始",
    "Srsf10": "精原细胞/起始", "Kpna2": "精原细胞/起始", "Utf1": "精原细胞/起始",
    "Tet1": "精原细胞/起始", "Dnmt3c": "精原细胞/起始", "Dnmt3l": "精原细胞/起始",
    "Msl3": "精原细胞/起始", "Pld6": "精原细胞/起始", "Boll": "精原细胞/起始",
    "Boule": "精原细胞/起始", "Morc1": "精原细胞/起始", "Eif4e": "精原细胞/起始",
    "Ddx20": "精原细胞/起始", "Sin3a": "精原细胞/起始",
    # 减数分裂前期
    "Spo11": "减数分裂前期", "Dmc1": "减数分裂前期", "Stra8": "减数分裂前期",
    "Sycp1": "减数分裂前期", "Sycp3": "减数分裂前期", "Rec8": "减数分裂前期",
    "Hormad1": "减数分裂前期", "Hormad2": "减数分裂前期", "Mei1": "减数分裂前期",
    "Mei4": "减数分裂前期", "Msh4": "减数分裂前期", "Msh5": "减数分裂前期",
    "Mlh1": "减数分裂前期", "Mlh3": "减数分裂前期", "Atm": "减数分裂前期",
    "Brca1": "减数分裂前期", "Brca2": "减数分裂前期", "Rad51": "减数分裂前期",
    "Rad50": "减数分裂前期", "Ankrd31": "减数分裂前期", "Fignl1": "减数分裂前期",
    "Tex11": "减数分裂前期", "Tex15": "减数分裂前期", "Meioc": "减数分裂前期",
    "Srsf1": "减数分裂前期", "Bcas2": "减数分裂前期", "Uba2": "减数分裂前期",
    "Ppp2ca": "减数分裂前期", "Baz1a": "减数分裂前期", "Upf2": "减数分裂前期",
    "Cnot4": "减数分裂前期", "Smg6": "减数分裂前期", "Mtr4": "减数分裂前期",
    "Hsf5": "减数分裂前期", "Iqch": "减数分裂前期", "Ptbp2": "减数分裂前期",
    # 精子形成（减数分裂后）
    "Ddx25": "精子形成", "Grh": "精子形成", "Miwi": "精子形成", "Piwil1": "精子形成",
    "Piwil4": "精子形成", "Mov10l1": "精子形成", "Tdrd": "精子形成", "Mael": "精子形成",
    "Brdt": "精子形成", "Crem": "精子形成", "Tnp1": "精子形成", "Tnp2": "精子形成",
    "Prm1": "精子形成", "Prm2": "精子形成", "Hspa2": "精子形成", "Oaz3": "精子形成",
    "Catsper": "精子形成", "Tssk": "精子形成", "Ift20": "精子形成", "Agfg1": "精子形成",
    # 体细胞
    "Arhgef15": "体细胞(Sertoli/Leydig)", "Nr2f2": "体细胞(Sertoli/Leydig)",
    "Lhcgr": "体细胞(Sertoli/Leydig)", "Insl3": "体细胞(Sertoli/Leydig)",
    "Gja1": "体细胞(Sertoli/Leydig)", "Nr5a2": "体细胞(Sertoli/Leydig)",
}

NODE_KW = [
    ("性别反转/转分化", re.compile(r"sex reversal|\bovarian\b|\bovary\b|transdifferentiat|\bSry\b|Foxl2|female germ", re.I)),
    ("环境/化学/物理", re.compile(r"busulfan|irradiat|heat stress|heat shock|cryptorchid|hyperthermia|cisplatin|cyclophosphamide|doxorubicin|phthalate|bisphenol|\bTCDD\b|dioxin|toxican|zinc deficien|folate deficien|vitamin A deficien|\bVAD\b|WIN ?18446", re.I)),
    ("激素/内分泌", re.compile(r"GnRH|hypophys|testosterone|androgen|\bFSH\b|gonadotropin|contracep|hormonal|estrogen", re.I)),
    ("精子形成", re.compile(r"spermiogenesis|spermatid|protamine|\bTnp[12]\b|\bPrm[12]\b|acrosome|flagell|chromatin remodel|elongat|round spermatid|residual body|cytoplasmic droplet", re.I)),
    ("减数分裂前期", re.compile(r"meiotic|meiosis|pachytene|zygotene|leptotene|diplotene|prophase|synapsis|recombination|crossover|double.?strand break|spermatocyte", re.I)),
    ("精原细胞/起始", re.compile(r"spermatogoni|stem cell|\bGDNF\b|\bPlzf\b|Zbtb16|Gfra1|Sall4|Dmrt1|Nanos[23]|Utf1|undifferentiated|self-renewal|progenitor", re.I)),
    ("体细胞(Sertoli/Leydig)", re.compile(r"sertoli|leydig|peritubular|blood-testis barrier|myoid", re.I)),
]

# 非转录组测定，要单独标出。规则要求明确的测定语，避免误伤
# （比如 "chromatin remodeler CHD8" 是转录组研究，不是 ChIP）
ASSAY_RULES = [
    ("ChIP-seq / 染色质结合", re.compile(r"\bChIP\b|ChIP-?seq|ChIPseq|CUT&RUN|CUT&Tag|\bATAC\b|H3K9me3|H3K4me3|H3K27ac|H3K36me3|histone modification|chromatin immunoprecipitation", re.I)),
    ("小 RNA 测序", re.compile(r"small ?RNA|microRNA (expression|profil|sequenc)|miRNA (expression|profil|sequenc)|piRNA (expression|profil|sequenc)|sRNA-?seq", re.I)),
    ("核糖体印迹 / 多聚核糖体", re.compile(r"ribosome profiling|ribo-?seq|RiboTag|polysom|translating ribosome", re.I)),
    ("甲基化 / 其他非表达组", re.compile(r"MeDIP|bisulfite|RRBS|WGBS|SNP array|comparative genomic hybrid|methylation profil", re.I)),
]

TISSUE_RULES = [
    ("分选生殖细胞", re.compile(r"GFP\+|GFP-|sorted|purified|FACS|flow.?sorted|isolated (germ|spermat|cell)|enriched", re.I)),
    ("单细胞", re.compile(r"single.?cell|scRNA|10x|drop.?seq|scRNA-seq", re.I)),
    ("特定细胞类型", re.compile(r"spermatocyte|spermatid|spermatogoni|sertoli|leydig|round spermatid|pachytene|leptotene", re.I)),
]


def node_of(text, gse_hint):
    for name, rx in NODE_KW:
        if rx.search(text):
            return name, "文中有阶段词"
    for g, n in GENE_NODE.items():
        if re.search(r"\b" + re.escape(g) + r"\b", text, re.I):
            return n, "按基因初判"
    return "待定", "需读文献"


def assay_of(text, kind):
    """先看有没有明确的其他测定语，没有就按 gdstype 判（比从标题猜可靠）"""
    for name, rx in ASSAY_RULES:
        if rx.search(text):
            return name
    if kind == "RNA-seq":
        return "转录组(RNA-seq)"
    if kind == "microarray":
        return "转录组(芯片)"
    return "待查"


def tissue_of(text):
    for name, rx in TISSUE_RULES:
        if rx.search(text):
            return name
    if re.search(r"whole testis|testes|testicular", text, re.I):
        return "全睾丸"
    return "待查"


def model_of(text):
    if re.search(r"knockout|\bKO\b|null|deficien|mutant|deletion|knock-down|knockdown|ablation|conditional", text, re.I):
        return "基因敲除/突变"
    if re.search(r"busulfan|cisplatin|cyclophosphamide|doxorubicin|phthalate|bisphenol|TCDD|dioxin", text, re.I):
        return "化学处理"
    if re.search(r"vitamin A|retinoic acid|WIN ?18446|zinc deficien|folate deficien|diet", text, re.I):
        return "营养/维A酸"
    if re.search(r"irradiat|heat|cryptorchid|hyperthermia", text, re.I):
        return "物理/环境"
    if re.search(r"GnRH|hypophys|testosterone|androgen|estrogen|hormonal", text, re.I):
        return "激素/内分泌"
    return "待定"


def raw_of(kind, suppfile, supp_state):
    sup = (suppfile or "").upper()
    if kind == "RNA-seq":
        return "SRA 有原始读长"
    if kind == "microarray":
        if re.search(r"CEL|IDAT|RAW|TAR|GPR", sup):
            return "有芯片原始文件"
        return "仅信号值"
    return "待查"


def matrix_of(supp_state):
    return {"matrix": "有（补充文件含矩阵）", "signal_only": "无（仅信号值）",
            "none": "无", "archive": "有（打包）", "other": "待查"}.get(supp_state, "待查")


def clean(s):
    """markdown 表格里不能有裸竖线和换行"""
    return str(s).replace("|", "／").replace("\n", " ").strip()


def main():
    meta = {r["gse"]: r for r in csv.DictReader(open(os.path.join(DATA, "geo_meta.tsv"), encoding="utf-8"), delimiter="\t")}
    cand = list(csv.DictReader(open(os.path.join(DATA, "geo_candidates.tsv"), encoding="utf-8"), delimiter="\t"))

    rows = []
    for c in cand:
        if c["bucket"] == "excluded":
            continue
        m = meta.get(c["gse"], {})
        text = " ".join([c.get("title", ""), m.get("summary", ""), c.get("sample_titles", ""), c.get("note", "")])
        node, basis = node_of(text, c["gse"])
        rows.append({
            "bucket": c["bucket"],
            "gse": c["gse"],
            "node": node,
            "node_basis": basis,
            "assay": assay_of(text, c["kind"]),
            "tissue": tissue_of(text),
            "model": model_of(text),
            "title": c["title"],
            "n_samples": c["n_samples"],
            "n_ko": c["n_ko_titles"],
            "n_wt": c["n_wt_titles"],
            "kind": c["kind"],
            "platform": ";".join("GPL" + p.strip() for p in str(c["platform"]).split(";") if p.strip()),
            "platformtitle": m.get("platformtitle", ""),
            "supp_state": c["supp_state"],
            "suppfile": c["suppfile"],
            "raw": raw_of(c["kind"], c["suppfile"], c["supp_state"]),
            "matrix": matrix_of(c["supp_state"]),
            "pdat": c["pdat"][:10],
            "pubmed": c["pubmed"],
            "sample_titles": c.get("sample_titles", "")[:200],
            "summary": m.get("summary", "")[:300],
        })

    cols = ["bucket", "gse", "node", "node_basis", "assay", "tissue", "model", "usability", "title",
            "n_samples", "n_ko", "n_wt", "kind", "platform", "platformtitle",
            "supp_state", "suppfile", "raw", "matrix", "pdat", "pubmed",
            "sample_titles", "summary"]
    for r in rows:
        if r["assay"].startswith("转录组"):
            r["usability"] = "可直接用" if r["matrix"].startswith("有") else "需从原始数据重算"
        else:
            r["usability"] = "非转录组，不纳入共性分析"
    with open(os.path.join(DATA, "inventory.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(rows)
    print("写入 data/inventory.csv（%d 条）" % len(rows))

    # ---------- 统计 ----------
    mainrows = [r for r in rows if r["bucket"] == "main"]
    print("\nmain %d 条，atlas %d 条" % (len(mainrows), len(rows) - len(mainrows)))
    print("\nmain 的测定类型:")
    for k, v in Counter(r["assay"] for r in mainrows).most_common():
        print("   %-16s %d" % (k, v))
    print("\nmain 的节点归类:")
    for k, v in Counter(r["node"] for r in mainrows).most_common():
        print("   %-22s %d" % (k, v))
    print("\nmain 的节点依据:")
    for k, v in Counter(r["node_basis"] for r in mainrows).most_common():
        print("   %-12s %d" % (k, v))
    print("\nmain 的组织类型:")
    for k, v in Counter(r["tissue"] for r in mainrows).most_common():
        print("   %-14s %d" % (k, v))
    print("\nmain 的原始数据可得性:")
    for k, v in Counter(r["raw"] for r in mainrows).most_common():
        print("   %-18s %d" % (k, v))

    # ---------- 生成 markdown ----------
    order = ["精原细胞/起始", "减数分裂前期", "精子形成", "体细胞(Sertoli/Leydig)",
             "激素/内分泌", "环境/化学/物理", "性别反转/转分化", "待定"]
    tx = [r for r in mainrows if r["assay"].startswith("转录组")]
    nontx = [r for r in mainrows if not r["assay"].startswith("转录组")]
    atlas = [r for r in rows if r["bucket"] == "atlas"]

    grp = defaultdict(list)
    for r in tx:
        grp[r["node"]].append(r)

    L = []
    L.append("# 小鼠生精阻滞模型 —— GEO 数据集清单\n")
    L.append("生成流程：`scripts/01_geo_search.py`（检索）→ `scripts/02_geo_meta.py`（取元数据、"
             "识别 KO/WT 样本、判断补充文件）→ `scripts/03_curate.py`（归类、出清单）。"
             "机读版见 `data/inventory.csv`。\n")

    L.append("## 一、总览\n")
    L.append("| 环节 | 数量 |")
    L.append("|---|---|")
    L.append("| GEO 检索命中（去重） | 673 |")
    L.append("| 其中小鼠数据集 | 449 |")
    L.append("| 睾丸相关 + KO/WT 样本各 ≥2（main） | %d |" % len(mainrows))
    L.append("| ├ **真转录组**（RNA-seq 或芯片） | **%d** |" % len(tx))
    L.append("| │  ├ 补充文件含表达矩阵，可直接用 | %d |" % len([r for r in tx if r["matrix"].startswith("有")]))
    L.append("| │  └ 无矩阵，需从 SRA 原始读长重算 | %d |" % len([r for r in tx if not r["matrix"].startswith("有")]))
    L.append("| └ 非转录组（ChIP / 小RNA / 核糖体 / 甲基化） | %d |" % len(nontx))
    L.append("| 发育图谱 / 细胞类型参考（另计） | %d |" % len(atlas))
    L.append("")

    L.append("## 二、可用的阻滞模型（只列转录组，按节点分组）\n")
    L.append("> **节点是推断的，不是 GEO 里的字段。** GEO 元数据只有标题、摘要、样本名，"
             "没有「阻滞在哪个阶段」。表中节点依据标题/摘要里的阶段词或基因名推断；"
             "标「按基因初判」的要读原文确认表型，标「需读文献」的连基因线索都没有。\n")

    for node in order:
        items = grp.get(node)
        if not items:
            continue
        L.append("\n### %s（%d 个）\n" % (node, len(items)))
        L.append("| GSE | 组织 | 模型 | 可用性 | 样本 | KO/WT | 平台 | 节点依据 | 年份 | 标题 |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for r in sorted(items, key=lambda x: x["pdat"]):
            L.append("| [%s](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s) | %s | %s | %s | %s | %s/%s | %s | %s | %s | %s |"
                     % (r["gse"], r["gse"], clean(clean(r["tissue"])), clean(clean(r["model"])), clean(clean(r["usability"])),
                        r["n_samples"], r["n_ko"], r["n_wt"], clean(clean(r["platform"])),
                        clean(clean(r["node_basis"])), r["pdat"][:4], clean(r["title"])[:58]))

    L.append("\n## 三、非转录组测定（%d 个，不能混进共性分析）\n" % len(nontx))
    L.append("这些看着像「睾丸 KO vs WT」，但测的不是表达量。如果按关键词一股脑抓下来"
             "直接和 RNA-seq 合并，结论会是错的。留在这里是为了说明筛掉了什么。\n")
    L.append("| GSE | 测定 | 节点 | 样本 | 平台 | 年份 | 标题 |")
    L.append("|---|---|---|---|---|---|---|")
    for r in sorted(nontx, key=lambda x: x["pdat"]):
        L.append("| %s | %s | %s | %s | %s | %s | %s |"
                 % (r["gse"], clean(clean(r["assay"])), clean(clean(r["node"])), r["n_samples"], clean(clean(r["platform"])),
                    r["pdat"][:4], clean(r["title"])[:58]))

    L.append("\n## 四、图谱 / 参考数据集（%d 个）\n" % len(atlas))
    L.append("发育图谱或细胞类型参考，用途是给 bulk 数据做去卷积，**不是阻滞模型本身**，"
             "不参与「有多少个独立模型」的计数。\n")
    L.append("| GSE | 测定 | 组织 | 样本 | 平台 | 年份 | 标题 |")
    L.append("|---|---|---|---|---|---|---|")
    for r in sorted(atlas, key=lambda x: x["pdat"]):
        L.append("| %s | %s | %s | %s | %s | %s | %s |"
                 % (r["gse"], clean(clean(r["assay"])), clean(clean(r["tissue"])), r["n_samples"], clean(clean(r["platform"])),
                    r["pdat"][:4], clean(r["title"])[:58]))

    L.append("\n---\n")
    L.append("## 五、这张清单说明了什么\n")
    L.append("1. **数据量够。** %d 个真转录组模型，三个节点都有覆盖，"
             "精原细胞/起始 %d 个、减数分裂前期 %d 个、精子形成 %d 个。"
             "原计划书里「凑不出 3 个 miRNA 敲除」的问题，换成阻滞模型框架就不存在了。\n"
             % (len(tx), len(grp.get("精原细胞/起始", [])),
                len(grp.get("减数分裂前期", [])), len(grp.get("精子形成", []))))
    L.append("2. **节点归类要靠读文献，不能靠元数据。** %d 个转录组模型里，"
             "只有 %d 个标题或摘要明确写了阶段；%d 个只能按基因名初判；%d 个连基因线索都没有。"
             "这是第一步真正的成本，不是技术问题。\n"
             % (len(tx), len([r for r in tx if r["node_basis"] == "文中有阶段词"]),
                len([r for r in tx if r["node_basis"] == "按基因初判"]),
                len([r for r in tx if r["node_basis"] == "需读文献"])))
    L.append("3. **关键词检索会混进非表达数据。** %d 个 main 里 %d 个不是转录组（占 %.0f%%），"
             "包括 ChIP-seq、小 RNA 测序、核糖体印迹、甲基化。这些在标题上都是「睾丸 + KO」，"
             "不筛就会混进来。\n" % (len(mainrows), len(nontx), 100.0 * len(nontx) / len(mainrows)))
    L.append("4. **组织类型不统一。** %s。不同组织（全睾丸 vs 分选细胞）的差异倍数不能直接比，"
             "这是后面必须处理的第一件事。\n"
             % "；".join("%s %d 个" % (k, v) for k, v in
                        Counter(r["tissue"] for r in tx).most_common()))
    L.append("5. **平台分散。** 模型跨 %d 个不同平台，批次效应是真实存在的，"
             "不是假想的风险。\n" % len({p for r in tx for p in r["platform"].split(";")}))

    with open(os.path.join(ROOT, "notes/数据集清单.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("\n写入 notes/数据集清单.md")
    print("  转录组 %d（可直接用 %d / 需重算 %d），非转录组 %d，图谱 %d"
          % (len(tx), len([r for r in tx if r["matrix"].startswith("有")]),
             len([r for r in tx if not r["matrix"].startswith("有")]), len(nontx), len(atlas)))


if __name__ == "__main__":
    main()
