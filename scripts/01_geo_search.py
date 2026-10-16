#!/usr/bin/env python3
"""
系统检索 GEO 上小鼠生精阻滞 / 雄性不育模型的数据集。

只做检索和元数据收集，不做任何分析。
输出 data/geo_hits.tsv 和 data/geo_detail.tsv。

用法：
    python3 scripts/01_geo_search.py
"""

import os
import json
import time
import urllib.request
import urllib.parse

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
os.makedirs(OUTDIR, exist_ok=True)

SLEEP = 0.4  # NCBI 无 key 限 3 req/s，留足余量


# ---------------------------------------------------------------- 检索词
# 分成几组，最后合并去重。每组都限定小鼠。
QUERIES = [
    # 核心表型词
    "spermatogenic arrest mouse testis",
    "meiotic arrest mouse testis",
    "arrest of spermatogenesis mouse",
    "male infertility mouse testis RNA-seq",
    "azoospermia mouse testis transcriptome",
    "oligozoospermia mouse testis",
    "germ cell depletion mouse testis",
    "Sertoli cell only mouse testis",
    "impaired spermatogenesis mouse RNA-seq",
    "defective spermatogenesis mouse testis",
    "spermatocyte arrest mouse",
    "round spermatid arrest mouse",

    # 减数分裂相关基因敲除
    "Spo11 knockout testis mouse",
    "Dmc1 knockout mouse testis",
    "Stra8 knockout mouse testis",
    "Sycp3 knockout mouse testis",
    "Sycp1 knockout mouse",
    "Rec8 knockout mouse testis",
    "Hormad1 knockout mouse",
    "Hormad2 knockout mouse testis",
    "Mei1 knockout mouse",
    "Mei4 knockout mouse testis",
    "Msh4 knockout mouse testis",
    "Msh5 knockout mouse",
    "Mlh1 knockout mouse testis",
    "Mlh3 knockout mouse testis",
    "Atm knockout mouse testis",
    "Brca1 knockout mouse testis",
    "Brca2 knockout mouse testis",
    "Rad51 knockout mouse testis",
    "Rad50 mouse germ cell",
    "Ankrd31 mouse testis",
    "Fignl1 mouse testis",
    "Tex11 knockout mouse testis",
    "Tex15 knockout mouse",
    "Hsf5 mouse testis",
    "Uba2 mouse testis spermatocyte",

    # 生殖细胞发育 / 精原细胞
    "Dmrt1 knockout mouse testis",
    "Dazl knockout mouse testis",
    "Boll knockout mouse testis",
    "Boule knockout mouse",
    "Nanos2 mouse testis",
    "Nanos3 mouse germ cell",
    "Sohlh1 knockout mouse testis",
    "Sohlh2 knockout mouse testis",
    "Plzf mouse spermatogonia",
    "Kit mouse spermatogonia",
    "Gfra1 mouse spermatogonia",
    "Neurog3 mouse testis",
    "Dcaf2 mouse spermatogonia",
    "Eif5 mouse germ cell",
    "20S proteasome mouse spermatogenesis",

    # piRNA / 小 RNA 通路
    "Piwil1 knockout mouse testis",
    "Piwil2 knockout mouse testis",
    "Mili mouse testis piRNA",
    "Miwi mouse testis",
    "Tdrd1 mouse testis",
    "Tdrd5 mouse testis",
    "Tdrd6 mouse testis",
    "Tdrd7 mouse testis",
    "Ddx4 mouse testis",
    "Mvh mouse testis",
    "Smg6 mouse testis piRNA",
    "Dnmt3l mouse testis",
    "Dnmt3a mouse germ cell",

    # miRNA 通路与具体 miRNA
    "microRNA knockout mouse testis",
    "miRNA knockout spermatogenesis mouse",
    "Dicer knockout mouse testis",
    "Dicer1 germ cell knockout mouse",
    "Dgcr8 knockout mouse testis",
    "Ago2 knockout mouse germ cell",
    "Drosha mouse germ cell",
    "miR-202 mouse testis",
    "miR-34b mouse testis",
    "miR-449 mouse testis",
    "miR-17-92 mouse gonad",
    "miR-506 mouse testis",
    "X-linked microRNA mouse testis",
    "miR-871 mouse testis",
    "miR-880 mouse testis",

    # 表观 / 染色质 / RNA 加工
    "H3mm mouse spermatogenesis",
    "histone variant mouse testis knockout",
    "chromatin remodeling mouse spermatocyte",
    "SRSF1 mouse male meiosis",
    "BCAS2 mouse testis",
    "Dhx36 mouse spermatogenesis",
    "Daxx mouse testis LINE1",
    "RNA polymerase II pausing mouse meiosis",
    "alternative splicing mouse spermatogenesis",
    "EWSR1 mouse testis meiosis",
    "CXXC1 mouse meiosis",
    "Msl3 mouse testis",
    "Setdb1 mouse germ cell",
    "Tet1 mouse germ cell",
    "Suv39h mouse testis",

    # 化学 / 环境 / 物理模型
    "busulfan mouse testis spermatogenesis",
    "busulfan mouse testis single cell",
    "cyclophosphamide mouse testis RNA-seq",
    "irradiation mouse testis spermatogenesis",
    "heat stress mouse testis spermatogenesis",
    "scrotal heat mouse testis",
    "varicocele mouse testis RNA-seq",
    "vitamin A deficiency mouse testis",
    "retinoic acid mouse testis spermatogenesis",
    "WIN 18446 mouse testis",
    "bisphenol mouse testis spermatogenesis",
    "phthalate mouse testis spermatogenesis",
    "DEHP mouse testis RNA-seq",
    "cadmium mouse testis RNA-seq",
    "triptolide mouse testis",
    "gossypol mouse testis",
    "doxorubicin mouse testis spermatogenesis",
    "cisplatin mouse testis spermatogenesis",
    "high fat diet mouse testis spermatogenesis",
    "diabetes mouse testis spermatogenesis RNA-seq",
    "obesity mouse testis spermatogenesis",
    "aging mouse testis spermatogenesis RNA-seq",
    "cryptorchidism mouse testis",

    # 单细胞参考图谱（用于去卷积）
    "mouse testis single cell atlas spermatogenesis",
    "mouse spermatogenesis single cell RNA-seq",
    "mouse testis single cell roadmap",
    "mouse spermatogonia single cell transcriptome",
    "mouse spermatocyte single cell RNA-seq",
]


def fetch(url, tries=4):
    """带退避的 GET。NCBI 偶尔会断连。"""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", "ignore")
        except Exception as e:
            if k == tries - 1:
                print("    放弃:", e)
                return None
            time.sleep(3 * (k + 1))
    return None


def esearch(term, retmax=30):
    url = (EUTILS + "/esearch.fcgi?db=gds&retmode=json&retmax=%d&term=%s"
           % (retmax, urllib.parse.quote(term)))
    raw = fetch(url)
    if not raw:
        return []
    try:
        return json.loads(raw)["esearchresult"]["idlist"]
    except Exception:
        return []


def esummary(ids):
    """esummary 一次最多 200 个 UID。"""
    out = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        url = (EUTILS + "/esummary.fcgi?db=gds&retmode=json&id=%s"
               % ",".join(chunk))
        raw = fetch(url)
        if not raw:
            continue
        try:
            res = json.loads(raw).get("result", {})
            for k, v in res.items():
                if isinstance(v, dict) and v.get("accession"):
                    out[k] = v
        except Exception as e:
            print("    esummary 解析失败:", e)
        time.sleep(SLEEP)
    return out


def main():
    print("=" * 60)
    print("第一步：跑检索词，收集 UID")
    print("=" * 60)

    uid2query = {}
    for i, q in enumerate(QUERIES, 1):
        ids = esearch(q)
        new = 0
        for uid in ids:
            if uid not in uid2query:
                uid2query[uid] = q
                new += 1
        print("  [%3d/%d] %-52s 命中 %2d 条 (新 %2d)"
              % (i, len(QUERIES), q[:52], len(ids), new))
        time.sleep(SLEEP)

    print("\n合计唯一 UID：%d 个" % len(uid2query))

    print("\n" + "=" * 60)
    print("第二步：拉元数据")
    print("=" * 60)
    detail = esummary(list(uid2query.keys()))
    print("拿到 %d 条记录" % len(detail))

    # ---------------- 写全量命中表 ----------------
    allrows = []
    for uid, d in detail.items():
        acc = d.get("accession", "")
        if not acc.startswith("GSE"):
            continue  # 跳过 GSM / GPL / GDS
        allrows.append({
            "gse": acc,
            "title": (d.get("title") or "").replace("\t", " ").replace("\n", " "),
            "gdstype": (d.get("gdstype") or "").replace("\t", " "),
            "taxon": d.get("taxon") or "",
            "n_samples": d.get("n_samples") or "",
            "platform": d.get("gpl") or "",
            "pdat": d.get("pdat") or "",
            "query": uid2query.get(uid, ""),
        })

    allrows.sort(key=lambda r: r["gse"])
    f_all = os.path.join(OUTDIR, "geo_hits.tsv")
    with open(f_all, "w", encoding="utf-8") as f:
        cols = ["gse", "title", "gdstype", "taxon", "n_samples",
                "platform", "pdat", "query"]
        f.write("\t".join(cols) + "\n")
        for r in allrows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("\n全量命中写入 %s（%d 条 GSE）" % (f_all, len(allrows)))

    # ---------------- 过滤：小鼠 + 睾丸相关 ----------------
    KEY_TISSUE = [
        "testis", "testes", "testicular", "spermatogen", "spermato",
        "germ cell", "germline", "germ line", "seminiferous",
        "spermatocyte", "spermatogoni", "spermatid", "meiotic",
        "meiosis", "male fertility", "male infertilit", "azoosperm",
        "oligozoosperm", "steril", "Sertoli", "Leydig",
    ]

    filtered = []
    for r in allrows:
        if r["taxon"] != "Mus musculus":
            continue
        blob = (r["title"] + " " + r["gdstype"]).lower()
        if not any(k.lower() in blob for k in KEY_TISSUE):
            continue
        filtered.append(r)

    f_filt = os.path.join(OUTDIR, "geo_detail.tsv")
    with open(f_filt, "w", encoding="utf-8") as f:
        cols = ["gse", "title", "gdstype", "n_samples", "platform", "pdat", "query"]
        f.write("\t".join(cols) + "\n")
        for r in filtered:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("过滤后（小鼠 + 睾丸相关）：%d 条，写入 %s" % (len(filtered), f_filt))
    print("\n完成。下一步跑 02_geo_supplementary.py 查原始数据可得性。")


if __name__ == "__main__":
    main()