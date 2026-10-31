# Targeted literature comparison and a paired-cytology candidate

## Scope and retrieval

On 2026-10-31, three explicit PubMed title/abstract queries targeted mouse arrest transcriptomes, comparative/integrative spermatogenesis transcriptomes, and testis cellular-composition/deconvolution studies. Their counts were 33, 224 and 45, respectively. All identifiers reconciled within each query, giving 296 unique PMID records with retrieved metadata and abstracts. Four directly relevant papers received focused abstract/method review; the remaining retrieved records were not all scientifically appraised. This is a targeted context search, **not a systematic review or proof of novelty**. Cross-study and composition queries deliberately admit other species, and the mouse-arrest query also retrieves some female-meiosis studies. Search matches are not eligibility decisions.

Queries, returned translations, ID sets, count reconciliation, response hashes, request failures and access dates are retained in local research provenance. All downloads used the configured proxy. The initial 100-record abstract batch first exceeded the response-size cap and then encountered TLS EOF; smaller 40-record batches completed through the same proxy. The adult TPM download also succeeded after a same-proxy TLS retry. There were 31 network attempts: 26 successes and five recorded failures. The candidate-acquisition follow-up exceeded its four-attempt local plan by one because the earlier retry was not counted when scheduling two metadata requests; this deviation is retained. No direct fallback was used. One requested full text returned HTTP500 twice and remains unavailable locally.

## Focused comparison

| Published work | Directly relevant approach | Content checked | Implication for this project |
|---|---|---|---|
| [Cappallo-Obermann et al., 2010](https://doi.org/10.1007/s00439-010-0865-9) | Human Y-chromosomal-deletion biopsies; reduce sample-composition differences and compare expression effects | DOI registration, PubMed abstract and primary full-text methods | Separating tissue-composition differences from expression effects is an established problem. Their sample selection is distinct from the present fixed-candidate gene matching. |
| [Margolin et al., 2014](https://doi.org/10.1186/1471-2164-15-39) | Mouse developmental whole-testis RNA-seq and computational decomposition of cell-associated expression | DOI registration, PubMed abstract and primary deconvolution methods | Developmental decomposition predates this project; its expression invariance and mixing assumptions require attention. It is not the present cross-perturbation partial-conjunction test. |
| [Ball et al., 2016](https://doi.org/10.1186/s12864-016-2865-1) | Matched cytological substage counts and RNA-seq from isolated juvenile germ-cell samples; permutation-based maximum covariance analysis | DOI registration, PubMed abstract, primary methods/results and GSE72833 source record | Paired cell-state measurements offer an orthogonal evaluation opportunity, with a restricted recovered-germ-cell denominator and source-specific exclusions. |
| [Hunnicutt et al., 2022](https://doi.org/10.1111/evo.14420) | Whole-testis versus enriched spermatogenic populations in house-mouse subspecies and fertile/sterile hybrids | DOI registration and PubMed abstract; primary full text not retrieved after two HTTP500 responses | The abstract directly reports composition-sensitive expression inference. Detailed design and quantitative benchmarks still require full-text review. |

These sources prevent presenting attention to composition, developmental staging or comparative expression as a unique contribution. The current contribution remains a scoped, reproducible analysis with explicit genetic-unit accounting, negative primary results, missingness and reference-transfer diagnostics. Whether that combination supports a publishable contribution remains an open review question. No assertion that prior studies lack these features is justified by four focused comparisons.

## GSE72833 candidate acquisition

The [GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72833) explicitly links Ball2016/PMID27519264 and lists 32 archive samples. The article reports 30 juvenile paired cytology/RNA samples, followed by omission of two P8 samples from subsequent analysis, leaving 28. It separately reports four adult pachytene/diplotene-enriched preparations, each pooled from six mice. These are source assertions, not a new certification of animals. The juvenile cells were isolated germ cells rather than unseparated whole testis.

Downloaded and hash-recorded:

- The source-labeled gene-expression TPM table: 20,368 rows and 28 expression columns, plus gene and substage fields.
- The adult STA-PUT TPM table: 85,537 feature rows and four columns; the gene identifier is an implicit leading field absent from the header labels.
- One juvenile GSM processing record, and the publisher's Tables S1-S10 archive. The archive was inventoried; its Excel tables have not yet been extracted or reconciled to GEO samples.

**Quantitative use is blocked pending unit clarification.** The juvenile table's observed values include entries above one million and column sums far above one million, inconsistent with ordinary nonnegative TPM in an unrepeated feature universe. The source explicitly calls this file RSEM TPM, so the conflict is recorded rather than silently rescaling or relabeling it. The adult table has a different numerical profile, close to one million per column. These observations do not establish a particular alternative unit or prove a specific author error.

The candidate is not added to the primary family or any expression benchmark. Next steps are to reconcile feature identity, deposited scale, the cytology table and native sample labels, preserve the authors' exclusion history, and then specify a held-out evaluation before inspecting its outcomes. Relative cell frequencies among recovered germ cells cannot certify absolute whole-testis fractions, and RNA contributions need not equal cell counts.

## Local audit files

`data/paired_cytology_candidate.json` records the acquisition/QC gate and input identities. Detailed search contracts, raw payloads and review notes remain ignored local research material. No new manuscript, raw data or literature attachment was uploaded, submitted or imported into Zotero in this task.
