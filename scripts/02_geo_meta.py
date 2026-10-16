#!/usr/bin/env python3
"""
对 01_geo_search.py 筛出的小鼠睾丸相关 GSE，逐条拉完整元数据。

关键字段（esummary 实际返回名）：
    accession, title, summary, gdstype, taxon, n_samples,
    gpl, platformtitle, pdat, suppfile, samples, pubmedids, ftplink, bioproject

输出：
    data/geo_meta.tsv        全量元数据
    data/geo_candidates.tsv  按用途分类后的候选清单

用法：
    python3 scripts/02_geo_meta.py
"""

import os
import re
import json
import time
import urllib.request
import urllib.parse

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SLEEP = 0.4
LIMIT = int(os.environ.get("GEO_LIMIT", "0"))


def fetch(url, tries=4, timeout=60):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-meta/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "ignore")
        except Exception as e:
            if k == tries - 1:
                print("    fail:", str(e)[:60])
                return None
            time.sleep(2 * (k + 1))
    return None


def gse_to_uid(gse):
    """GSE63166 -> 200063166；GSE126936 -> 200126936。数字补足 6 位。"""
    num = gse[3:]
    return "200" + num.zfill(6)


def esummary_batch(uids):
    url = (EUTILS + "/esummary.fcgi?db=gds&retmode=json&id=" + ",".join(uids))
    raw = fetch(url)
    if not raw:
        return {}
    try:
        res = json.loads(raw).get("result", {})
    except Exception as e:
        print("    parse fail:", e)
        return {}
    out = {}
    for k, v in res.items():
        if isinstance(v, dict) and v.get("accession"):
            out[v["accession"]] = v
    return out


# ---------------------------------------------------------------- 判定逻辑
SEQ_HINT = ["high throughput sequencing", "rna-seq", "rna seq"]
ARRAY_HINT = ["array", "microarray"]
OTHER_HINT = ["genome binding", "methylation", "non-coding rna", "third-party", "other"]

MATRIX_SUFFIX = {"TXT", "TSV", "CSV", "MTX", "XLSX", "XLS", "COUNT", "COUNTS", "FPKM", "TPM"}
SIGNAL_SUFFIX = {"BW", "BED", "WIG", "BIGWIG"}          # 信号文件，不是计数矩阵
ARCHIVE_SUFFIX = {"TAR", "GZ", "ZIP", "RDS", "H5"}       # 需解开或专用读取


def classify_kind(gdstype, platformtitle):
    blob = ((gdstype or "") + " " + (platformtitle or "")).lower()
    if any(h in blob for h in SEQ_HINT):
        return "RNA-seq"
    if any(h in blob for h in ARRAY_HINT):
        return "microarray"
    if "genome binding" in blob:
        return "ChIP/ATAC"
    if "methylation" in blob:
        return "methylation"
    if "non-coding" in blob:
        return "smallRNA"
    return "other"


def classify_supp(suppfile):
    """判断补充文件里有没有可能直接用的表达矩阵。"""
    if not suppfile:
        return "none", ""
    parts = [p.strip().upper() for p in re.split(r"[,;\s]+", suppfile) if p.strip()]
    has_matrix = bool(set(parts) & MATRIX_SUFFIX)
    has_signal = bool(set(parts) & SIGNAL_SUFFIX)
    has_arch = bool(set(parts) & ARCHIVE_SUFFIX)
    if has_matrix:
        return "matrix", ",".join(parts)
    if has_arch:
        return "archive", ",".join(parts)
    if has_signal:
        return "signal_only", ",".join(parts)
    return "other", ",".join(parts)


KO_PAT = re.compile(
    r"(\bko\b|knockout|knock-out|\bcKO\b|null|mutant|\b-/-\b|delet|\bcKO\b|"
    r"\bKO[-_ ]|conditional|ablat|deplet)", re.I)
WT_PAT = re.compile(r"(\bwt\b|wild[- ]?type|control|ctrl|\bhet\b|flox|\blox\b|sham|vehicle|mock)", re.I)


def group_signal(titles):
    """从样本标题里找 KO 和对照的线索。"""
    ko = [t for t in titles if KO_PAT.search(t)]
    wt = [t for t in titles if WT_PAT.search(t)]
    return len(ko), len(wt), "; ".join(titles[:8])


TISSUE_PAT = re.compile(
    r"(testis|testes|testicular|spermat|germ ?cell|germline|seminiferous|"
    r"spermatocyte|spermatogoni|spermatid|meiotic|meiosis|Sertoli|Leydig|azoosperm)", re.I)
ATLAS_PAT = re.compile(r"(atlas|roadmap|single[- ]cell|scRNA|trajector|develop|landscape|map of)", re.I)


# ---------------------------------------------------------------- 主流程
def main():
    src = os.path.join(DATA, "geo_detail.tsv")
    if not os.path.exists(src):
        print("找不到", src, "，先跑 01_geo_search.py")
        return

    rows = []
    with open(src, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    gses = [r["gse"] for r in rows if r.get("gse")]
    if LIMIT:
        gses = gses[:LIMIT]
    print("待拉元数据 %d 条" % len(gses))

    meta = {}
    B = 50
    for i in range(0, len(gses), B):
        chunk = gses[i:i + B]
        uids = [gse_to_uid(g) for g in chunk]
        got = esummary_batch(uids)
        meta.update(got)
        print("  [%3d/%3d] 本批拿到 %d" % (min(i + B, len(gses)), len(gses), len(got)), flush=True)
        time.sleep(SLEEP)

    # ---------------- 整理 ----------------
    full = []
    for gse in gses:
        d = meta.get(gse)
        if not d:
            continue
        titles = [s.get("title", "") for s in (d.get("samples") or [])]
        nko, nwt, sample_blob = group_signal(titles)
        supp_state, supp_detail = classify_supp(d.get("suppfile", ""))
        kind = classify_kind(d.get("gdstype", ""), d.get("platformtitle", ""))
        blob = (d.get("title", "") + " " + (d.get("summary") or "")).lower()

        full.append({
            "gse": gse,
            "title": (d.get("title") or "").replace("\t", " ").replace("\n", " ")[:200],
            "kind": kind,
            "platform": d.get("gpl", ""),
            "platformtitle": (d.get("platformtitle") or "")[:80],
            "n_samples": d.get("n_samples", ""),
            "suppfile": supp_detail,
            "supp_state": supp_state,
            "n_ko_titles": nko,
            "n_wt_titles": nwt,
            "pdat": d.get("pdat", ""),
            "pubmed": ",".join(d.get("pubmedids") or []),
            "tissue_hit": "yes" if TISSUE_PAT.search(blob) else "no",
            "atlas_hint": "yes" if ATLAS_PAT.search(blob) else "no",
            "sample_titles": sample_blob[:280],
            "summary": (d.get("summary") or "").replace("\t", " ")[:300],
        })

    f_full = os.path.join(DATA, "geo_meta.tsv")
    cols = ["gse", "title", "kind", "platform", "platformtitle", "n_samples",
            "suppfile", "supp_state", "n_ko_titles", "n_wt_titles", "pdat",
            "pubmed", "tissue_hit", "atlas_hint", "sample_titles", "summary"]
    with open(f_full, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for rec in full:
            f.write("\t".join(str(rec[c]).replace("\t", " ") for c in cols) + "\n")
    print("\n全量写入 %s（%d 条）" % (f_full, len(full)))

    # ---------------- 分类 ----------------
    main_cand, atlas_cand, excluded = [], [], []
    for rec in full:
        if rec["tissue_hit"] != "yes":
            excluded.append((rec, "非睾丸组织"))
            continue
        if rec["kind"] == "RNA-seq" and rec["n_ko_titles"] >= 2 and rec["n_wt_titles"] >= 2:
            if rec["atlas_hint"] == "yes" and rec["n_ko_titles"] < 3:
                atlas_cand.append((rec, "图谱/发育"))
            else:
                main_cand.append((rec, ""))
        elif rec["kind"] == "RNA-seq" and rec["atlas_hint"] == "yes":
            atlas_cand.append((rec, "图谱"))
        elif rec["kind"] == "microarray" and rec["n_ko_titles"] >= 2 and rec["n_wt_titles"] >= 2:
            main_cand.append((rec, "芯片，作验证"))
        else:
            excluded.append((rec, "类型/分组不明确"))

    f_cand = os.path.join(DATA, "geo_candidates.tsv")
    with open(f_cand, "w", encoding="utf-8") as f:
        f.write("\t".join(["bucket", "gse", "title", "kind", "platform", "n_samples",
                           "supp_state", "suppfile", "n_ko_titles", "n_wt_titles",
                           "pdat", "pubmed", "note", "sample_titles"]) + "\n")
        for bucket, items in [("main", main_cand), ("atlas", atlas_cand), ("excluded", excluded)]:
            for rec, note in items:
                f.write("\t".join(str(x).replace("\t", " ") for x in [
                    bucket, rec["gse"], rec["title"], rec["kind"], rec["platform"],
                    rec["n_samples"], rec["supp_state"], rec["suppfile"],
                    rec["n_ko_titles"], rec["n_wt_titles"], rec["pdat"],
                    rec["pubmed"], note, rec["sample_titles"]]) + "\n")

    print("\n分类结果：")
    print("  主分析候选 %d" % len(main_cand))
    print("  参考图谱候选 %d" % len(atlas_cand))
    print("  排除 %d" % len(excluded))
    print("\n写入 %s" % f_cand)

    # 顺手回答几个已知问题
    print("\n" + "=" * 50)
    print("平台核查（之前问过的几个）")
    for g in ["GSE3436", "GSE23119", "GSE112393", "GSE126936"]:
        for rec in full:
            if rec["gse"] == g:
                print("  %-11s %-12s %-26s n=%-5s supp=%-8s %s"
                      % (g, rec["kind"], rec["platformtitle"][:26],
                         rec["n_samples"], rec["supp_state"], rec["title"][:50]))


if __name__ == "__main__":
    main()
