# Cross-model analysis of spermatogenic-arrest transcriptomes: historical local execution plan

> **Historical design, superseded by the workflow audited on 2026-10-22. This is not the current execution plan.** References below to a passed positive control, cell fractions, percentages attributable to composition, or residual “true signal” are retained as historical statements, not validated conclusions. See the [current entry point](README.md), [methods and limitations](docs/METHODS.md), and [audited results](docs/RESULTS.md). This English translation preserves the earlier design and its numerical records so that its errors and evolution remain traceable. Historical timings, environment measurements, and proposed acceptance criteria are not new verification results.

> Historical target machine: Apple M4 / 10 cores / **16 GB RAM** / approximately 80 GB free disk / macOS aarch64. This document originally replaced `notes/原计划书.md`. Local feasibility records were documented in `docs/01-environment.md` through `docs/04-pipeline-runbook.md`.

## 1. Original rationale

The early workspace contained two conflicting documents:

| Document | Historical assessment | Reason |
|---|---|---|
| `notes/原计划书.md` (214 lines) | Unsupported by the available data | Identifying core genes shared by at least three miRNA knockouts requires at least three corresponding transcriptomes. The historical search found only miR-202; n=1 cannot support cross-model intersections or similarity matrices. |
| `notes/思路整理.md` (219 lines) | Adopted research direction | The question changed from finding shared signals to asking whether those signals are reproducible and how much they could reflect changes in cellular composition. The original claim that its biological foundation was validated was subsequently qualified. |

The plan adopted the research question from the second document and reported an end-to-end positive-control run on the GSE236388 XAP5 cKO testis arm: fitted germ-cell weight decreased by **0.194**, fitted somatic weight increased by **0.194**, and reconstruction R² was **0.680–0.742**, against a threshold of 0.6. The original plan described this as a passed positive control and evidence of a biologically expected direction. **These observations demonstrate workflow operation and directional consistency, not validated cell-number fractions.** The historical command uses `code/06_validate_deconvolution.py`; original output is recorded in `docs/02-datasets.md`.

The proposed design used eight whole-testis bulk datasets plus the GSE112393 single-cell reference: first compare cross-model log fold changes, then construct a composition-only null model, and finally examine residual expression. The original machine-time estimate was about one hour; human work was estimated at **8–15 hours** for metadata checks, figure review, and reporting.

## 2. Why the original miRNA proposal was rejected

| Original assumption | Historical audit finding |
|---|---|
| Six classes of miRNA-knockout mouse testis transcriptomes were available | Only miR-202 was identified: GSE174731 / GSE126936 / GSE99160. |
| miR-34b/c+miR-449, miR-506, and XmiRs datasets were available | Matching datasets were not found in that search; this is not a claim of universal absence. |
| miR-17-92 datasets GSE225675/GSE225676 were suitable | These represent XY complete sex reversal and were considered unsuitable comparisons for spermatogenic arrest. |
| A Jaccard similarity matrix and clustering could be constructed | A single perturbation cannot support a cross-model intersection or similarity matrix. |
| References 3 and 5 were distinct | Both listed the same Wu et al. 2014 PLoS Genetics paper. |

Two additional concerns motivated the revised question:

1. Mixing RNA-seq and microarrays can introduce platform effects that resemble cross-model agreement.
2. Bulk-testis consensus patterns, including lower `Prm1/Prm2/Tnp1/Tnp2` and higher `Sox9/Wt1/Gma4/Amh/Cyp11a1`, may reflect the loss of later germ-cell populations. The historical plan treated this as a composition hypothesis to test. It does not establish that composition fully explains these genes or every arrest model; the original gene spelling is retained.

## 3. The proposed three-part design

```text
Part 1  Describe cross-model convergence in log fold changes.
        Reproduce the observation that some arrest models appear similar.
Part 2  Construct a composition-only null model.
        Ask how much cross-model similarity the fitted composition predicts.
Part 3  Examine residual patterns after removing the prediction.
        Residual patterns are not automatically cell-intrinsic biological signal.
```

### Mathematical definition of the historical null model

The predicted expression contrast propagates a composition shift through a reference profile and applies the logarithm:

$$\hat c_m(g)=\log\!\big((S\,\bar p^{KO}_m)(g)+1\big)-\log\!\big((S\,\bar p^{WT}_m)(g)+1\big).$$

It asks which genes would appear changed in model m if only composition changed while each cell population's expression profile remained fixed. The reference's applicability and the fitted weights remain assumptions.

The historical simulations reported three implementation concerns:

1. Projecting logFC directly onto the column space of S returned **0.585** when the target was **1.0**, an underestimate of about **40%**, even with a matched reference. The plan attributed this to the nonlinearity of the logarithm; using the group-average weights propagated through S returned **0.992** in that simulation.
2. A naive Mantel test returned the same **p = 0.158** under both a composition-only artifact scenario and a scenario containing additional signal. The plan replaced it with multiple regression on distance matrices (MRM), using model pairs. Its historical rationale cited Legendre 2015, the mismatch between sample-distance matrices and within-model KO-versus-WT contrasts, and the simulation's lack of discrimination. This rationale is not a newly verified literature claim.
3. Permutation was defined over **models**, not genes or samples, by pairing model i's logFC with model j's predicted composition contrast. Gene shuffling produced a narrow null standard deviation of **0.007**, compared with **0.129** under model reassignment in the historical simulation. Proposed exhaustive counts were **720** for M=6, **5040** for M=7, and preferably **40320** for M=8.

The proposed main regression was:

$$y_{ij}=\mathrm{cov}(\mathrm{lfc}_i,\mathrm{lfc}_j)\sim \mathrm{cov}(\hat c_i,\hat c_j)+\mathbb 1[\text{node}_i=\text{node}_j]+\cdots.$$

### Power limitations recorded in the plan

| M (number of models) | a=0, composition-only artifact | a=0.20 | a=0.40 |
|---|---|---|---|
| 6 | Unstable: 1.59, **sd 19.9** | 0.293 ± 0.158 | 0.103 ± 0.061 |
| 12 | — | sd decreased to 0.03 | — |

The plan required reporting compShare as a point estimate with a **95% permutation interval**, rather than a seemingly precise two-decimal percentage. It argued that M=6–8 could address whether composition is a major explanation, but could not resolve **37% versus 44%**. This was a historical design aspiration, not a guarantee of adequate power or valid composition attribution; the limitation must remain explicit.

### Scientific wording and prior work

The proposed framing was “test the assumptions underlying this type of analysis,” rather than claiming that earlier authors were wrong.

The plan identified **Zhang J et al., “Microarray meta-analysis identifies candidate genes for human spermatogenic arrest,” Andrologia 2019;51(8):e13301**, PMID **31037746**, DOI **10.1111/and.13301**. Its historical notes reported a Crossref title check and described 178 human testes, seven microarray datasets, and four categories based on Johnsen scores. These are retained source notes, not a new literature audit.

The full text was not obtained because it was behind the Wiley paywall; the historical record reported no PMC copy and Europe PMC `isOpenAccess:N`. The abstract did not discuss composition/deconvolution, but this does not establish what the inaccessible full text addressed. **Do not claim that this paper ignored composition.** The recorded first author was Zhang J, whereas `notes/思路整理.md` used “Kui F.”

## 4. Historical CORE-8 dataset inventory

The historical workflow restricted deconvolution to samples labeled `whole_testis`; `code/05_deconvolution.R` stops for other tissue classes. Sorted populations were excluded from this composition interpretation because their mixture is determined by the sorting design. Their fitted weights were not treated as in vivo biological composition estimates.

| # | GSE | Model | Historical node | Total n (KO/WT) | Age | Platform | Historical processing status |
|---|---|---|---|---|---|---|---|
| 1 | GSE38096 | Germ-cell Sin3a cKO, Vasa-cre | Spermatogonial/onset | 6 (2/2 + 2 HET) | P0 | GPL13912 Agilent 4x44K | Normalization required |
| 2 | GSE42601 | Dazap1-Fn, arrest before MI | Spermatogonial/onset | 6 (3/3) | P21 | GPL1261 Mouse430_2 | RMA required |
| 3 | GSE294209 | Rnu11 cKO | Meiotic prophase | 6 (3/3) | P14 | GPL24247 NovaSeq | TPM, not counts |
| 4 | GSE236388 | XAP5 cKO, D16 testis arm only | Meiotic prophase | 6 (3/3) | D16 | GPL24247 NovaSeq | Preferred count input |
| 5 | GSE101511 | Hnrnpgt cKO | Meiotic prophase | 6 (3/3) | 18 dpp | GPL19057 NextSeq | Transcript-to-gene mapping required |
| 6 | GSE283904 | N6amt1-SKO | Meiotic prophase | 8 (4/4) | Adult in the historical table; currently unconfirmed | GPL24247 NovaSeq | Preferred count input |
| 7 | GSE68283 | Rfx2-/- | Postmeiotic/spermiogenic | 9 (4/5) | P29–30 | GPL17021 HiSeq | Preferred count input |
| 8 | GSE68284 | Rfx2-/-, P21 arm of the same paper | Postmeiotic/spermiogenic | 11 (6/5) | P21 | GPL17021 HiSeq | Preferred count input |

Historical node coverage was **2 / 4 / 2**, and **8/8** cohorts were labeled `tissue_class=whole_testis`. The plan listed five deconvolution inputs as “yes” and three as “yes-with-caveat,” with zero “no.” This historical eligibility label does not establish current quantitative validity. Sample-source statements are retained in `data/dataset_selection.csv`, column `whole_testis_evidence`. The current baseline excludes the two Sin3a HET samples, does not accept an unverified adult label for N6amt1, and does not count the two Rfx2 ages as independent genetic models.

### Inventory corrections recorded in the plan

The plan reported **15 mislabeled entries** in `data/inventory.csv` and warned against using that inventory without checking source records.

Entries labeled whole testis that did not represent that material included:

- **GSE323391:** testis organoids in vitro.
- **GSE237030:** 10X single-cell RNA data and ES-induced somatic cells.
- **GSE125495:** female gonadal somatic cells.
- **GSE310340:** although the series theme concerned JQ1 inhibition of meiotic prophase, the 20 GSMs under this accession were recorded as **Kidney×10 + Liver×10**, from `Sample_source_name_ch1`. The plan regarded this as a likely toxicity-arm subseries and excluded it as a testis source.
- **GSE23119:** sorted spermatogonia, BSA gradient, purity >90%.
- **GSE236388:** n=12 combined six D16 testes with six D56 cauda epididymal sperm profiles; only the six testis profiles were eligible.
- **GSE221433:** whole testis, but a Western-diet comparison rather than a genetic arrest model; assigned a control role.

Other recorded corrections:

- GSE158471 / GSE253154 / GSE288152 / GSE324232 were sorted-cell sources despite unresolved inventory labels.
- `supp_state=matrix` was incorrect for GSE259417 and GSE114042, which provided author DEG tables, and GSE57745, which provided transposable-element-family counts rather than gene counts.
- GSE235777 contained **3KO/3WT**, not the inventory's 3/6; GSE101511 contained **3/3**, absent from the original inventory.
- **GSE259417** provided only an author DESeq2 DEG table, also missing `KO_rep2`. The plan allowed examination of provided logFC but rejected deconvolution; it did not establish an independent expression matrix.

Historical EXTENDED-4 reserve sources were **GSE213442 (Kpna2)**, **GSE214201 (Iqch, matrix 4WT/3KO)**, **GSE235777 (Msl3)**, and **GSE221433 (Tas1r3, dietary control)**.

### Historical download budget

| Input scope | Bytes | At 126 KB/s |
|---|---|---|
| CORE-8 supplementary inputs | 121,550,497 | 0.26 h |
| GSE112393 reference | 82,756,793 | 0.18 h |
| Total | **204,307,290** | Approximately **27 minutes** |

The plan excluded `GSE244681_RAW.tar`, recorded as **462 GB**, which would require about **44 days** at 126 KB/s. These are historical estimates, not current transfer measurements.

## 5. Historical local feasibility and resource constraints

### 5.1 Dense-reference memory risk

The recorded GSE112393 matrix dimensions were **37,241 genes × 34,633 cells**, obtained from the compressed-file header.

| Approach | Recorded memory | Historical assessment |
|---|---|---|
| Dense float64 | 10.32 GB; copies during transpose, scaling, or subsetting could raise use to 20.6 GB | Rejected |
| Dense float32 | 5.16 GB | Little headroom |
| Chunks of 1000 gene rows | 139 MB per chunk | Feasible but unnecessary for aggregation |
| Stream directly into reference S | Measured peak **47.7 MB** | Adopted |

The recorded failure example was:

```r
G <- 37241; N <- 34633
m <- matrix(0, nrow=G, ncol=N)   # Succeeded; object.size reported 9.6 Gb.
m2 <- m * 2                    # cannot allocate vector of size 9.6 Gb
                               # vector memory limit of 16.0 Gb reached, see mem.maxVSize()
```

The historical R vector-heap limit was `mem.maxVSize()` = **16,384 Mb = 16 GiB**, distinct from the operating system's memory policy. Streaming avoids the demonstrated allocation/copy problem; the example does not establish that all memory failures have the same cause.

Other recorded sizes were: 5,000 reference genes, **1.39 GB**; 2,000 genes, **554 MB**; streamed S of **37,241×11**, **3.3 MB**; one bulk matrix of **20,000×40**, **6.4 MB**; measured full DESeq2 peak, **0.93 GB**. The plan prohibited eight simultaneous DESeq2 jobs, budgeting **8 × 2 GB = 16 GB**.

### 5.2 Historical environment choice

- The system R library was not writable (`Operation not permitted`), so the plan used a project-local `.Rlib` through `.libPaths(c(".Rlib", .libPaths()))`. It recorded **171 packages / 934 MB**, with `library()` passing for **21** key packages.
- CRAN and Bioconductor **3.23** supplied arm64 binaries from `bin/macosx/sonoma-arm64/contrib/4.6/`. The 21 checked packages showed no compilation traces and did not require clang in that installation; this corrected an earlier claim that precompiled binaries were unavailable.
- Historical installation attempts for `MuSiC` / `bisque` failed after six GitHub retries. The plan used NNLS with `nnls` and `quadprog`, with installed `EpiDISH` as an alternative.
- A scanpy wheel was available for Python **3.14.7**. The local `.venv` contained numpy, pandas, scipy, sklearn, openpyxl, and matplotlib.
- Conda was **not universally unusable**. The historical A/B environment check succeeded after redirecting `XDG_CACHE_HOME` and `CONDA_PKGS_DIRS` and setting `CONDA_REGISTER_ENVS=false`. The earlier failure concerned writing `~/.conda/environments.txt`. The plan retained venv + `.Rlib` as an engineering choice; the workaround is in `docs/01-environment.md`.
- Direct PyPI access had produced connection resets in the historical environment; `envs/setup.sh` used the Tsinghua mirror, with a recorded rate of about **2.6 MB/s**. This does not replace the current local-proxy download requirement.
- Matplotlib used `MPLCONFIGDIR="$PWD/.mplcache"` because `~/.matplotlib` was not writable.

### 5.3 Historical download implementation

Recorded GEO FTP throughput was **15–35 KB/s** for a single stream. Aggregate measurements were: one stream **21 KB/s**, four **61 KB/s**, eight **126 KB/s**, and sixteen **81 KB/s**. Eight streams were preferred for that measurement; this is not a current concurrency recommendation for every endpoint.

`aria2c`, `wget`, `axel`, and `lftp` were unavailable. The plan used `scripts/fetch_geo.py`, a standard-library range downloader with eight concurrent requests, **3 MB** chunks, and chunk-level resume.

Two failures and fixes were recorded:

1. Large ranges generated by macOS/BSD `seq` could be formatted in scientific notation, causing `value too great for base (error token is "2e")` and zero-byte files. Chunk lists were generated in Python instead.
2. Holding chunks in memory until assembly lost a 20-minute download after a transient `HTTP 503`. The downloader was changed to save each chunk immediately, retry with exponential backoff, and call `rmtree(.parts)` after successful assembly. Stale incomplete parts had overwritten a complete reference, changing **81,261,540 B** to **79,322,404 B**. `gzip_ok()` added full CRC/ISIZE verification and skipped already complete inputs.

## 6. Historical stepwise execution proposal

The plan estimated **40–60 minutes** of machine time with a global peak below **1.5 GB**. These were design estimates, not a guarantee for the current workflow. Each proposed step wrote outputs and was intended to be independently resumable.

| Step | Proposed task | Outputs | Estimated machine time | Resume |
|---|---|---|---|---|
| 0 | Environment checks: .Rlib/.venv/disk/memory/network, 11 assertions | Check report | <1 min | — |
| 1 | Download CORE-8 and reference; verify size, gzip, SHA-256 | `data/raw/`, `_manifest.tsv` | ~30 min, background | Yes |
| 2 | Per-dataset input and gene-space harmonization: normalize SYMBOLs, filter low expression, exclude `^mt-`; prohibit cross-GSE cbind/merge with stop assertions | `data/interim/<GSE>_expr.rds` | ~5 min | Yes |
| 3 | Separate DE per dataset: counts with limma-voom or DESeq2, microarrays with limma; fixed columns `gene, logFC, AveExpr, t, P.Value, adj.P.Val` | `results/de/<GSE>_de.csv` | ~10 min | Yes |
| 4 | Part 1: cross-model logFC similarity matrix, heatmap, clustering | `results/act1/` | ~2 min | Yes |
| 5 | Stream S, fit P with NNLS, and assess 80/20 held-out-gene R² | `data/ref/`, `results/deconvolution/` | ~12 min plus 727 s to construct S | Yes |
| 6 | Part 2: MRM and model-reassignment permutations | `results/act2/` | ~5 min | Yes |
| 7 | Part 3: residualization, residual DE, and cross-model comparison; use limma rather than DESeq2 because residuals can be negative | `results/act3/` | ~10 min | Yes |
| 8 | Calibration and sensitivity using known-mixture pseudobulk controls | `results/calibration/` | ~10 min | Yes |
| 9 | Figures fig01–fig15 and tables table01–table10 | `results/figures/`, `results/tables/` | ~5 min | Yes |

Historical human-time estimates were **1–2 h** for metadata/group checks, **4–8 h** for figure review and reporting, and **8–15 h** overall. Suggested scheduling was seven sessions at two hours per week, approximately **3.5 weeks**, around the author's examination/interview schedule. These personal scheduling notes are not part of the current scientific execution contract.

### Historical Step 5 acceptance criteria

- S should have **37,241 × 11** entries, with genes in rows; every column of P should sum to 1.
- Median reconstruction R² should exceed **0.6**; the recorded GSE236388 values were **0.680–0.742**.
- The KO contrast should shift fitted germ-cell weight downward and somatic weight upward; recorded changes were **−0.194 / +0.194**. The original “POSITIVE CONTROL PASS” label concerns that workflow check, not absolute fraction validity.
- S should be CPM-normalized. Historical simulation notes reported somatic-weight overestimation of **0.26** without normalization, in the same direction as the composition-shift hypothesis. CPM reduced the recorded error from **0.151** to **0.033**, a **4.6×** improvement.
- Recorded raw mean-count row sums by class were: meiotic **10,576**, postmeiotic **6,058**, somatic **2,418**, and spermatogonial **2,179**; meiotic was **4.9** times spermatogonial. These observations motivated normalization, not proof of transfer to real tissues.

### Three historical reference-parsing pitfalls

1. The first `DGE.txt.gz` row contained **34,633** cell barcodes without a gene-name column; data rows contained **34,634** fields. The recorded pandas settings were `header=None, skiprows=1, index_col=0`. A scalar dtype caused `could not convert string to float: '0610005C13Rik'`; the workaround used `{i: np.float32 for i in range(1, ncell+1)}`.
2. `PerCellAttributes` had CRLF endings, a header beginning with `#`, and no final newline. `wc -l` returned **34636**; `awk 'NR>5'` returned **34632**, one below the actual **34633** cell labels. Use parsed `len(labels)` rather than line-count assumptions.
3. The series mixed whole-testis and enriched/sorted material. ST1–6 had enrichment `none`; other entries included diploid/interstitial/sertoli/spermatogonia. Filtering using `PerCellAttributes` was required for the intended source scope.

## 7. Uncertainties recorded in the historical plan

| ID | Concern | Proposed response |
|---|---|---|
| U-1 | Low power: at M=6–8, a true compShare near 0.75 may be difficult to distinguish from 0.90 | Use MRM rather than relying on an unstable ratio denominator; report intervals and explicit limitations. |
| U-2 | Adult mouse reference versus juvenile bulk, including GSE38096 P0 and GSE68284 P21, or different strains | Consider systematic transfer bias, report age differences, and assess age-stratified sensitivities. |
| U-3 | Only 161/220 node labels had explicit literature statements; restricting to reliable nodes reduces M further | Make the reliability/power tradeoff explicit. |
| U-5 | Circular calibration: generating and deconvolving pseudobulk with GSE112393 does not test reference compatibility with real tissue | Add cross-paradigm controls, such as the sorted resource PMC4837615. The historical expectation that these controls might fail was a hypothesis, not a result. |
| U-7 | The proposed concern that NNLS orthogonality mechanically lowers residual similarity was rejected in the recorded simulation | Residual similarity instead increased at a=0: same-node versus cross-node **δ=+0.092, z=+2.57**. The plan proposed subtracting this artifact baseline; this is not an established correction for all real datasets. |

### Historical calibration checklist

The original heading referred to **eight** criteria fixed before inspecting real data, but the listed items were **C1, C2, C3, C5, C6, C7, and C8**; no C4 result is supplied here. Preserve that gap rather than invent an eighth result.

| Criterion | Historical recorded result |
|---|---|
| C1: mismatched-reference proportion-recovery MAE <0.05 | 0.024 |
| C2: compShare in [0.8,1.2] at a=0 | 1.020 |
| C3: monotonic decline with added signal | 1.02→0.82→0.50→0.21→0.07 |
| C5: orthogonality does not mechanically suppress residual similarity | Difference −0.005 |
| C6: CPM improves on raw profiles by at least 2× | 4.6× |
| C7: permutation null is nondegenerate | Reported as passing |
| C8: held-out R² is discriminative | Reported as passing |

These simulation checks do not establish biological calibration or transferability.

## 8. Decisions proposed at the time

1. **Node reliability versus power:** only **161/220** labels had explicit supporting literature. The plan proposed reliable nodes for the primary analysis and all nodes for a sensitivity, recognizing that M=6–8 was already a power limitation.
2. **Unknown cells:** the recorded Unknown fraction was **6.37%**. Excluding these cells produced a recovery error of **0.0496**, versus **0.0205** when assigned to Postmeiotic, approximately **2.4** times worse, with bias aligned to the hypothesis. The historical proposal was an exclusion main branch plus an explicit inclusion sensitivity, rather than silent reassignment; that proposal is not validation of the reassignment.
3. **Cal-1 circularity:** same-source pseudobulk calibration cannot establish external compatibility. Cross-paradigm controls were required, with conclusion strength reduced if they failed.

## 9. Historical risk ranking

| Priority | Risk | Proposed mitigation |
|---|---|---|
| 1 | Dense Step 5 allocation and `vector memory limit of 16.0 Gb reached` | Stream reference aggregation; recorded peak 47.7 MB. |
| 2 | Reversing WT/KO labels may not raise an error but reverses the contrast | Step 2/3 assertions and 1–2 h of manual source checks. |
| 3 | `GSE244681_RAW.tar`, 462 GB / 44 days at the recorded rate | Exclude that asset from the local download plan. |
| 4 | Eight concurrent DESeq2 jobs, 8×2 GB = 16 GB | Run sequentially. |
| 5 | NCBI HTTP 503 / throttling | Resumable chunks, retries, exponential backoff. |
| 6 | Wrong matrix orientation or log scale | Explicit orientation checks; the historical notes suggested inspecting gene-name rows with `grep`. |
| 7 | GSE259417 author DEG table without expression matrix | Restrict to supplied logFC examination; no deconvolution. |

## 10. Historical file map

Original source filenames are retained so that existing paths remain usable; Chinese text within filenames is not a separate language edition.

```text
spermatogenic-arrest-plan/
├── README.md                         Current repository entry point, English by default
├── README.zh-CN.md                   Chinese README with reciprocal language link
├── PLAN.md                           This historical plan, translated into English
├── LICENSE
├── .gitignore
├── docs/
│   ├── 01-environment.md              Historical 335-line environment report
│   ├── 02-datasets.md                 Historical 230-line dataset audit and budget
│   ├── 03-null-model.md               Historical 843-line null-model specification
│   └── 04-pipeline-runbook.md         Historical 913-line execution/failure runbook
├── code/
│   ├── 04_normalize_signature.py     S to CPM; row-sum/column-name checks
│   ├── 05_deconvolution.R            NNLS, tissue_class gate, held-out R²
│   └── 06_validate_deconvolution.py  Historical end-to-end positive-control checker
├── scripts/
│   ├── 01_geo_search.py … 03_curate.py  Inventory generation
│   ├── build_signature.py             Streaming S construction, 47.7 MB peak
│   └── fetch_geo.py                   Resumable range download and gzip checks
├── envs/                              setup.sh / install_pkgs.R / _drive.sh
├── data/
│   ├── README.md                      Original-data retrieval instructions
│   ├── inventory.csv                  Historical 220-row inventory; 15 reported errors
│   ├── dataset_selection.csv          Historical 57 rows × 19 columns:
│   │                                  core 8 / extended 4 / rejected 30 / mislabel 15
│   └── samples.tsv                    Historical 86 sample rows, including tissue_class
└── notes/
    ├── 原计划书.md                     Original rejected proposal
    ├── 思路整理.md                     Source of the revised research question
    └── 数据集清单.md                   Historical dataset list
```

Runtime paths excluded from version control included `.Rlib/` (historically 171 packages / 934 MB), `.venv/`, `.cache/`, `logs/`, `data/raw/`, `data/ref/`, `data/interim/`, `envs/vendor/`, and `envs/probe/`. Original data are not included in the repository.

## 11. Historical first-run example

Use the current [README](README.md) for supported execution. The following command is preserved for reproducing the earlier positive-control diagnostic, not for validating cell fractions or running every current analysis branch.

```bash
cd spermatogenic-arrest-plan

# 0) Supply the two required files; original data are not committed.
#    See data/README.md.
#    data/ref/signature_GSE112393_cpm.csv.gz
#    data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz

# 1) Historical idempotent environment check, approximately 15 s.
bash envs/setup.sh

# 2) Reproduce the earlier positive-control diagnostic, approximately 1 min.
MPLCONFIGDIR="$PWD/.mplcache" .venv/bin/python code/06_validate_deconvolution.py \
    --sig data/ref/signature_GSE112393_cpm.csv.gz \
    --bulk data/raw/GSE236388/GSE236388_XAP5_counts.txt.gz \
    --wt-pattern WT --ko-pattern cKO
# Historical output: POSITIVE CONTROL PASS; germ-cell delta = -0.194;
# R² 0.680–0.742. This output is not biological calibration.
```
