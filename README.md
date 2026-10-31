**English** | [简体中文](README.zh-CN.md)

# Cross-model transcriptomic analysis of spermatogenic arrest

This project compares shared and heterogeneous expression patterns in mouse whole-testis transcriptomes and tests whether single-cell references support interpretations based on cellular composition. It reanalyzes public data; no new animal or sequencing experiments were performed.

The frozen baseline is an exploratory analysis with audited data and methods: **8 cohorts, 56 main-analysis samples, and 7 genetic perturbations**. Agreement with an expected direction is not treated as validation of cell-fraction accuracy, and residual expression is not equated with a shared molecular mechanism.

## Main findings

- Some model pairs show stable positive expression-effect correlations, whereas others show clear inverse relationships. The evidence does not support treating every arrest model as the same transcriptomic response.
- All 56 individual sample omissions, together with TPM-pseudocount and probe-summary alternatives, produced 63 sensitivity variants. Several stronger associations retained their direction; the weaker Hnrnpgt associations were unstable. Samples were not excluded on this basis.
- Single-cell references mismatch early meiotic populations in independent sorted-cell data. Resolving the reference into 12 germ-cell states did not adequately resolve this problem. Validated cell-number fractions and causal percentages attributable to composition are therefore not reported.

![KO-minus-WT effect correlations across genes tested in all eight cohorts.](results_snapshot/all_cohort_comparison/eight_cohort_effects.png)

The figure shows Spearman correlations over 8,099 commonly tested genes, not significance tests based on independent biological replicates. TPM and microarray analyses are sensitivity branches and are not fully interchangeable with count-based effects. See [results and limitations](docs/RESULTS.md) and the [compact result snapshot](results_snapshot/README.md).

## Optional data expansion (2026-10-23)

The original eight-cohort baseline remains frozen. A targeted exploratory branch reviewed 32 GEO series and incorporated Dcaf2 P9, Pramef12 P2/P7 with heterozygous controls, and Msl3 P45 as a phenotype-negative comparison. This added 24 whole-testis samples and separately examined 20 sorted WT expression profiles. These distinct roles are not pooled as additional arrest cases. Methods, file-identity issues, and reproduction scope are documented in the [data expansion report](docs/DATA_EXPANSION.md); outputs are in `expansion_snapshot/`.

## Subsequent analyses and frozen baseline

**2026-10-26 — [Setx developmental and Eef1g sorted-cell source audit](docs/NEW_MODEL_CONTEXTS.md).** Two GEO queries returned 61 distinct series; 12 sources were reviewed. The branch downloaded and separately analyzed 40 Setx whole-testis TPM profiles and 6 Eef1g sorted L/Z count profiles. Among the fixed 36 downward priorities, 24, 19, and 33 retained a same-direction change of at least 1.5-fold under all individual sample omissions at Setx P14, P16, and P18, respectively. Cellular composition remains a possible explanation; these sources add no votes to the primary test. Fresh-directory reproduction matched 14 tables, 2 JSON records, and 2 figures. Outputs: `new_models_snapshot/`; entry point: `scripts/run_new_models.py`.

**2026-10-26 — [Directional replication across genetic models and candidate applicability](docs/REPLICABILITY.md).** Five count cohorts containing 40 samples were tested as four genetic units; the two Rfx2 ages were not counted twice. The complete directional hypothesis family over 9,675 genes produced zero discoveries under both BY and BH adjustment. The 105 descriptive candidates and 36 downward priorities stable under individual sample omissions are reported separately, without being described as confirmed shared mechanisms. Stage disagreements between two single-cell resources and post-hoc coverage checks in six existing additional contrasts are retained. Fresh-directory reproduction matched 28 tables, 4 JSON records, and 3 figures. Outputs: `replicability_snapshot/`; entry point: `scripts/run_replicability.py`.

**2026-10-25 — [Leptotene reference limits and independent UMI reconciliation](docs/REFERENCE_LIMITS.md).** The original GSE107644 matrix and author tables S1/S2 were reconciled for 1,136 cells with explicit stage labels, with additional gene-set, weighting, length, and source-library holdout sensitivities. Adding L/Z reference profiles did not reliably resolve external mismatches and could assign spermatogonial controls mainly to the meiotic class. Rechecking the original twenty sorted profiles did not resolve Cnot4 LZ. Outputs remain separate in `reference_limits_snapshot/`; the candidate reference is not promoted to a calibrated cell-fraction model. Check: `scripts/check_reference_limits.py`.

**2026-10-25 — [Independent leptotene-source diagnostics](docs/LEPTOTENE_FOLLOWUP.md).** Four candidate studies were reviewed, incorporating six stage-separated RNA count profiles from GSE155967 and four mixed LZ/pachytene controls from GSE83264. The two explicitly labeled leptotene profiles were still assigned mainly to spermatogonia by the finer reference; the two zygotene profiles shifted to the broad meiotic class. Synchronization, age, sorting, genotype, and library-technical limitations are retained. PreL is excluded from broad-label scoring. No new reference was trained from these results, and the eight-cohort baseline was not changed. Outputs: `leptotene_snapshot/`; check: `scripts/check_leptotene.py`.

**2026-10-25 — [Complete P15 EmptyDrops early-cell diagnostics](docs/EMPTYDROPS_P15.md).** All 13,168 author-labeled cells were reconciled. The finer reference's improvement extended to all 938 zygotene cells, including 672 newly incorporated cells, but 907 leptotene cells were still assigned mainly to spermatogonia. Early meiotic mismatch is therefore not considered resolved. Aggregation and fitting agreed between the two count versions when barcodes, labels, and genes were matched. Fresh-directory reproduction matched 14 tables and two design/provenance JSON records. Outputs remain separate in `emptydrops_p15_snapshot/`; check: `scripts/check_emptydrops_p15.py`.

**2026-10-24 — [Reference compatibility with explicit author labels in three libraries](docs/LABELED_REFERENCE.md).** The adult B6, P5, and P15 libraries were reconciled for 8,791 cells, with diagnostics for 38 author-labeled type groups and identification of explicit leptotene/zygotene labels. An independent sensitivity examined 266 overlapping zygotene cells with cached expression. Expression for most low-complexity cells had not yet been retrieved in this branch, so it did not establish resolution of LZ mismatch or cell-fraction accuracy. Outputs: `labeled_reference_snapshot/`; check: `scripts/check_labeled_snapshot.py`.

**2026-10-24 — [Developmental-source audit and fixed-reference transfer diagnostics](docs/REFERENCE_TRANSFER.md).** Eight capture libraries and 15,882 cells were audited. Held-out-gene prediction was evaluated over 20 existing sorted profiles and 32 library-processing combinations. The finer reference improved prediction error, without establishing cell-fraction accuracy; cell-level author labels were still missing in this branch. General offline integrity check: `scripts/check_project.py`.

**2026-10-24 — [Early Lrh1 contrasts and Ppp2ca file audit](docs/FOLLOWUP_LRH1.md).** Sixteen Lrh1 samples were separated into two comparisons, retaining author-batch adjustment and method sensitivities. Ppp2ca transcript FPKM was audited but not mixed into the count-based primary analysis. These exploratory results are separate in `followup_snapshot/`.

### Frozen eight-cohort baseline

| Cohort | Genetic perturbation | KO / WT | Data type |
|---|---|---:|---|
| GSE236388 | XAP5, D16 testis arm | 3 / 3 | Gene counts |
| GSE283904 | N6amt1 | 4 / 4 | Gene counts |
| GSE68283 | Rfx2, P29–30 | 4 / 5 | Gene counts |
| GSE68284 | Rfx2, P21 | 6 / 5 | Gene counts |
| GSE101511 | Hnrnpgt, 18 dpp | 3 / 3 | Legacy transcript-feature counts |
| GSE294209 | Rnu11, P14 | 3 / 3 | Author-provided TPM |
| GSE38096 | Sin3a, P0 | 2 / 2 | Agilent, Cy3 channel only |
| GSE42601 | Dazap1, P21 | 3 / 3 | Affymetrix mouse4302 |

The two Rfx2 age cohorts come from the same study and are not independent genetic models. The two Sin3a HET samples are excluded from the main contrast; the XAP5L D56 cauda epididymal sperm arm is excluded from the testis analysis. Current official sample records do not confirm the N6amt1 age, so the unverified historical adult label is not used for age analyses. Original histological evidence has not been established for every cohort's arrest stage; historical node labels do not support causal stage assignments.

Cohorts are processed independently. Gene-level effects, rather than raw expression matrices, are compared across cohorts. Composition exploration is restricted to the five count cohorts. Official URLs, sizes, and SHA-256 hashes for the 22 inputs are recorded in the [asset manifest](data/analysis_assets.json); source and group definitions are in the [data documentation](data/README.md).

## Reproduction

Python, R, and curl must already be installed. The recorded local runtime is Python 3.14.7 / R 4.6.1; different R/Bioconductor annotation versions may change gene mappings. Dependency installation does not require moving system interpreters or using sudo.

```bash
git clone https://github.com/Gensouno/spermatogenic-arrest-plan.git
cd spermatogenic-arrest-plan

# Initial installation; add --check to inspect an existing environment.
bash envs/setup_analysis.sh

# Download missing inputs and verify frozen checksums.
bash envs/python.sh scripts/data_assets.py --download

# Recompute the 19-stage baseline workflow from its inputs.
bash envs/python.sh scripts/run_analysis.py

# Verify the compact result snapshot included in the repository.
bash envs/python.sh scripts/export_snapshot.py --verify
```

Use `scripts/data_assets.py --verify` to check existing inputs, or `--import-from /path/to/existing-project` to import verified inputs from another directory. A checksum mismatch stops processing instead of silently overwriting an existing file.

### Optional download proxy

Copy `envs/download-proxy.example` to `.download-proxy.local` in the project root and replace its URL with the local HTTP(S) proxy address. This file is Git-ignored; do not commit proxy addresses or credentials. Downloads and dependency installation launched through `envs/python.sh` or `envs/rscript.sh` inherit the configuration. If the proxy fails, repair it instead of falling back silently to a direct connection. Without this file, the launchers preserve the caller's proxy environment. Direct system Python/R commands do not read this project-specific configuration.

A fresh local directory was also used to rerun the baseline without earlier outputs or intermediate files, and key result tables agreed. This reused installed dependencies: it is not a clean installation or a different-machine test, and it did not test every first-time download after deleting all original inputs. Timing and scope are documented in the [reproduction report](docs/REPRODUCIBILITY.md).


Check every available frozen branch (including local follow-ups) with `bash envs/python.sh scripts/check_all_snapshots.py`. This verifies artifact and bound-source integrity; numerical reproduction and biological validation are separate.

## Project navigation

- [Visual research progress](docs/RESEARCH_PROGRESS.md): completed computational work, current findings, manuscript status and the remaining evidence needed.

The root README defaults to English; use the language links above for the Chinese edition. Some linked technical reports and historical source filenames remain in Chinese.

- [Setx age-interaction analysis](docs/SETX_INTERACTION.md): joint age/genotype model, source-code and barcode limitations, sensitivities, and reproduction.
- [Setx source and length audit](docs/SETX_SOURCE_AND_LENGTH.md): archive identities, allele context, and length-aware count sensitivity.
- [Setx omission stability](docs/SETX_STABILITY.md): count omissions, matched-universe diagnostics, and explicit model dependence.
- [Additional whole-testis contexts](docs/INDEPENDENT_CONTEXTS.md): source audits and separate Rif1, Hormad1/2, Mael and Morc1 count analyses; the primary family remains frozen.
- [Reference-informed candidate matching](docs/STAGE_SPECIFICITY.md): abundance/stage-matched backgrounds, distinct stage descriptors, and the refreshed six-sample Rif1 age audit.
- [N6amt1 primary-source audit](docs/N6AMT1_SOURCE.md): eight reconciled source records, unresolved sampling age, and the distinction between CTL labels and verified wild-type genotypes.
- [Competitive priority-set testing](docs/COMPETITIVE_SET.md): all four additional contexts have a downward set direction, but none passes joint adjustment with estimated residual correlation; all pass under the assumed-correlation sensitivity. Results are model dependent and leave the original primary family unchanged.
- [Targeted literature comparison](docs/LITERATURE_CONTEXT.md): prior composition-aware studies and a newly acquired paired-cytology candidate; the candidate's deposited units require clarification before quantitative use.
- [English manuscript synchronization workflow](docs/PAPER_WORKFLOW.md): incremental writing, evidence registration, checks, Word export, and private-release boundaries.
- [Methods](docs/METHODS.md): platform-specific processing, gene universes, sensitivities, and composition diagnostics.
- [Results](docs/RESULTS.md): numerical findings, negative results, and interpretive limitations.
- [Reproducibility](docs/REPRODUCIBILITY.md): execution order, environment, and verification scope.
- [Interview discussion guide](docs/INTERVIEW.md): datasets, code, and limitations.
- [History](docs/HISTORY.md): earlier AI drafts, rejected assumptions, and reasons for retaining them.

Current implementations are in `code/` and `scripts/`. The baseline's compact, checksummed outputs are in `results_snapshot/`; complete regenerated outputs in `results/` are Git-ignored. Original data, dependency environments, and large working matrices are not committed.

`PLAN.md`, numbered reports `docs/01`–`04`, `notes/`, and the three legacy numbered analysis scripts are historical material, not the current execution plan. [PLAN.md](PLAN.md) is available in English and retains its historical status.

## Sources and licensing

The single-cell reference is [Green et al., Developmental Cell, 2018](https://pubmed.ncbi.nlm.nih.gov/30146481/); the external sorted-cell resource is [da Cruz et al., BMC Genomics, 2016](https://link.springer.com/article/10.1186/s12864-016-2618-1). Cohort-specific sources are recorded in the asset manifests and sample-level GEO evidence.

Original code is released under the [MIT License](LICENSE). Upstream datasets, papers, and supplementary materials retain their original licenses and citation requirements; this repository's software license does not replace them.
