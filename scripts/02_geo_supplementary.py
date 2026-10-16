#!/usr/bin/env python3
"""
对 01_geo_search.py 筛出来的 GSE，逐条查：

  1. 平台是芯片还是测序（决定要不要纳入主分析）
  2. 有没有原始计数矩阵 / supplementary 文件
  3. 原始 fastq 在不在 SRA（决定能不能重定量）
  4. 样本标题，用来判断分组和阻滞节点

输出 data/geo_probe.tsv。

用法：
    python3 scripts/02_geo_supplementary.py
"""

import os
import re
import json
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTPBASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
SLEEP = 0.4
LIMIT = int(os.environ.get("GEO_PROBE_LIMIT", "0"))  # 0 = 全跑


def fetch(url, tries=3, timeout=45):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-probe/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "ignore")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if k == tries - 1:
                return None
            time.sleep(2 * (k + 1))
        except Exception:
            if k == tries - 1:
                return None
            time.sleep(2 * (k + 1))
    return None


def probe_ftp(gse):
    """列出 GEO FTP 上这个 GSE 目录里的文件，判断数据可得性。"""
    stub = gse[:-3] + "nnn" if len(gse) > 6 else gse
    url = "%s/%s/%s/suppl/" % (FTPBASE, stub, gse)
    html = fetch(url)
    if html is None:
        # suppl 目录不存在，退回看 GSE 根目录
        url = "%s/%s/%s/" % (FTPBASE, stub, gse)
        html = fetch(url)
        if html is None:
            return {"suppl_files": "", "has_matrix": "no", "has_raw_counts": "unknown",
                    "has_soft": "no", "file_kinds": ""}
        names = re.findall(r'href="([^"?/][^"]*)"', html)
        return {"suppl_files": ";".join(names[:12]),
                "has_matrix": "no",
                "has_raw_counts": "unknown",
                "has_soft": "yes" if any(n.endswith(".soft.gz") for n in names) else "no",
                "file_kinds": ""}

    names = re.findall(r'href="([^"?/][^"]*)"', html)
    names = [n for n in names if not n.startswith("/")]

    kinds = set()
    for n in names:
        low = n.lower()
        if low.endswith((".txt.gz", ".tsv.gz", ".csv.gz", ".txt", ".tsv", ".csv")):
            kinds.add("table")
        elif low.endswith((".tar", ".tar.gz", ".tgz")):
            kinds.add("tarball")
        elif low.endswith((".xlsx", ".xls")):
            kinds.add("excel")
        elif low.endswith((".mtx.gz", ".mtx")):
            kinds.add("mtx")
        elif low.endswith(".h5"):
            kinds.add("h5")
        elif low.endswith((".bam", ".sam")):
            kinds.add("bam")
        elif low.endswith(".rds"):
            kinds.add("rds")

    joined = " ".join(names).lower()
    # 计数矩阵的常见命名
    count_hint = any(k in joined for k in [
        "count", "counts", "matrix", "umi", "expression", "fpkm", "rpkm",
        "tpm", "readcount", "gene_expression",
    ])
    # 只有差异分析结果、没有计数矩阵的情况
    deg_hint = any(k in joined for k in ["deg", "deseq", "edger", "diff", "vsd"])

    return {
        "suppl_files": ";".join(names[:12]),
        "has_matrix": "yes" if count_hint else ("deg_only" if deg_hint else "no"),
        "has_raw_counts": "likely" if ("table" in kinds and count_hint) else "unknown",
        "has_soft": "yes",
        "file_kinds": ",".join(sorted(kinds)),
    }


def probe_samples(gse):
    """用 GEO 的 series matrix 头拿样本标题和平台。"""
    stub = gse[:-3] + "nnn" if len(gse) > 6 else gse
    url = ("%s/%s/%s/matrix/%s_series_matrix.txt.gz"
           % (FTPBASE, stub, gse, gse))
    # 不下载 gz，改用 GEO 网页上的样本列表
    html = fetch("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief"
                 % gse, timeout=40)
    if not html:
        return {"platform": "", "n_samples": "", "sample_titles": "", "summary": ""}

    plat = re.findall(r"!Series_platform_id\s*=\s*(.+)", html)
    n = re.findall(r"!Series_sample_id\s*=\s*(.+)", html)
    summ = re.findall(r"!Series_summary\s*=\s*(.+)", html)
    design = re.findall(r"!Series_overall_design\s*=\s*(.+)", html)

    # 样本标题
    titles = []
    for line in html.split("\n"):
        m = re.match(r"!Sample_title\s*=\s*(.+)", line)
        if m:
            titles = [t.strip().strip('"') for t in m.group(1).split("\t")]
            break

    return {
        "platform": plat[0].strip() if plat else "",
        "n_samples": str(len(n)),
        "sample_titles": " | ".join(titles[:16]),
        "summary": (" ".join(summ + design))[:400].replace("\t", " "),
    }


PLATFORM_KIND = {
    "GPL1261": "microarray",
    "GPL6246": "microarray",
    "GPL8321": "microarray",
    "GPL6885": "microarray",
    "GPL6887": "microarray",
    "GPL16570": "microarray",
    "GPL21103": "RNA-seq",
    "GPL21273": "RNA-seq",
    "GPL24247": "RNA-seq",
    "GPL17021": "RNA-seq",
    "GPL19057": "RNA-seq",
    "GPL13112": "RNA-seq",
    "GPL28471": "RNA-seq",
    "GPL24676": "RNA-seq",
}


def main():
    src = os.path.join(DATA, "geo_detail.tsv")
    if not os.path.exists(src):
        print("找不到 %s，先跑 01_geo_search.py" % src)
        return

    rows = []
    with open(src, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    print("待核查 %d 条" % len(rows))
    if LIMIT:
        rows = rows[:LIMIT]
        print("（限制只跑前 %d 条）" % LIMIT)

    out = []
    for i, r in enumerate(rows, 1):
        gse = r["gse"]
        print("[%3d/%d] %s" % (i, len(rows), gse), end="  ", flush=True)

        ftp = probe_ftp(gse)
        time.sleep(SLEEP)
        smp = probe_samples(gse)
        time.sleep(SLEEP)

        plat = smp.get("platform", "") or r.get("platform", "")
        kind = PLATFORM_KIND.get(plat, "")
        if not kind:
            gdstype = (r.get("gdstype", "") or "").lower()
            if "high throughput sequencing" in gdstype:
                kind = "RNA-seq"
            elif "array" in gdstype:
                kind = "microarray"

        rec = {
            "gse": gse,
            "title": r["title"][:150],
            "kind": kind,
            "platform": plat,
            "n_samples": smp.get("n_samples") or r.get("n_samples", ""),
            "has_matrix": ftp["has_matrix"],
            "file_kinds": ftp["file_kinds"],
            "suppl_files": ftp["suppl_files"],
            "sample_titles": smp.get("sample_titles", "")[:300],
            "pdat": r.get("pdat", ""),
        }
        out.append(rec)
        print("%-12s n=%-4s matrix=%-8s files=%s"
              % (kind or "?", rec["n_samples"], rec["has_matrix"], rec["file_kinds"]))

    dst = os.path.join(DATA, "geo_probe.tsv")
    cols = ["gse", "title", "kind", "platform", "n_samples", "has_matrix",
            "file_kinds", "suppl_files", "sample_titles", "pdat"]
    with open(dst, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for rec in out:
            f.write("\t".join(str(rec[c]).replace("\t", " ") for c in cols) + "\n")

    print("\n写入 %s（%d 条）" % (dst, len(out)))

    # 简单汇总
    from collections import Counter
    print("\n平台分布:", dict(Counter(r["kind"] or "unknown" for r in out)))
    print("有计数矩阵:", dict(Counter(r["has_matrix"] for r in out)))


if __name__ == "__main__":
    main()