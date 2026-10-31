# N6amt1 primary-cohort source audit

The current [GSE283904 records](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283904)
were reviewed on 2026-10-31. Eight explicit GEO-to-BioSample/experiment relations were
reconciled through a count-first SRA lookup. The source contains four CTL and four
N6amt1-SKO testis libraries. All eight current GEO and SRA records omit sampling age.
The older local `adult` entries therefore remain unsupported curation, rather than
source-confirmed ages. Original input files and frozen results were preserved.

`CTL` is the literal source control-genotype label. The older `WT` analysis group
denotes that contrast's control arm; it does not establish a genetically wild-type
control. Exact floxed/Cre alleles, individual animals, pooling, batches and matched
histology have not been independently established. No expression model or primary
hypothesis family was changed by this audit.

A related published study was identified: Luo and colleagues, *Loss of N-6
adenine-specific DNA methyltransferase 1 leads to meiotic prophase abnormalities
and male sub-fertility in mice*, Biology of Reproduction (2025),
[DOI 10.1093/biolre/ioaf052](https://doi.org/10.1093/biolre/ioaf052),
[PMID40079145](https://pubmed.ncbi.nlm.nih.gov/40079145/). Its bibliographic identity
and abstract were checked. It reports meiotic abnormalities and reduced male fertility
following germ-cell N6amt1 deletion. This is related biological context: an explicit
GSE283904 accession link could not be verified from the accessible article content.
The publisher request returned403 and full methods were not retrieved. Neither the
article's animal ages nor an unverified dataset linkage was imputed into the cohort.

The initial proxy request failed with an SSL EOF; a single same-proxy retry succeeded.
Public source payloads, hashes and retrieval/access limitations are recorded in
`data/n6amt1_source_inputs.json`. Curated audit outputs exclude provider/contact fields.

```bash
# Requires the frozen cached payloads listed in data/n6amt1_source_inputs.json.
bash envs/python.sh scripts/audit_n6amt1_source.py
```

A fresh output directory reproduced the sample table and summary from the same cached
inputs. That establishes parsing and identity consistency, not animal independence,
phenotype validation or resolution of the missing sampling age. The primary zero
discovery result remains unchanged. The installed CAMERA function was inspected only
for technical availability; no competitive gene-set test was executed in this task.
