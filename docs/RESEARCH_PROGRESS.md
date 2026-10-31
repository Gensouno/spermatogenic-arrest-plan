# Research progress at 2026-10-31

The project has produced reproducible computational results, an English manuscript draft maintained locally, and versioned analysis code. A shared cell-intrinsic mechanism and submission readiness have not been established. Stage status is reported directly, without an estimated completion percentage.

![Research stages and their evidence boundaries](figures/progress-20261031/research_progress.png)

## What has been obtained

- **Audited comparisons and explicit negative primary results.** Five primary count cohorts contain 40 libraries representing four genetic units. The joint family of 38,700 directional hypotheses over 9,675 genes has zero BY/BH discoveries. The two Rfx2 ages contribute one unit. This limits a formal shared-expression claim; it does not demonstrate that shared biology is absent.
- **A focused descriptive candidate set.** The fixed 36 downward priorities are useful for follow-up, with measurement and omission-stability denominators retained. They are not a new FDR-controlled discovery set.
- **Evidence of source and model dependence.** Four additional count contexts retain different numbers of priorities. Reference-informed backgrounds narrow descriptive contrasts. Separate CAMERA tests are downward in all contexts, but adjusted support depends on estimated versus assumed residual correlation; test degrees of freedom also change.
- **Source constraints rather than imputed facts.** Refreshed Rif1 source records annotate P14; N6amt1 sampling age remains absent, and its CTL label is not proof of literal wild type. Setx source and animal/technical uncertainties remain explicit.
- **A better-defined literature and data backlog.** A targeted PubMed search retrieved 296 unique records; four papers received focused appraisal. Four new verified bibliographic identities bring the local bibliography to 20 entries. This is not a systematic review or proof of novelty. GSE72833 expression files and supplements were acquired, but juvenile scale/feature identity and cytology/sample mapping must be resolved before quantitative evaluation.

![Descriptive retention and competitive-test model dependence](figures/progress-20261031/research_findings.png)

## Practical use

The code and frozen artifacts let a collaborator reproduce the specified processed-input calculations and inspect missingness, source assumptions and model sensitivity. The descriptive priorities focus future analysis or experimental planning without identifying a proven mechanism. The English draft preserves the methods, results and limitations for eventual scientific review. Submission requires additional source verification, contribution assessment, authorship, ethics/licenses and venue-specific preparation.

The release checks recorded 142 passing tests, 18 consistent frozen snapshots, 336 frozen artifacts and 357 bound source/specification checks. The private manuscript check has zero errors and 36 draft/human-review/submission warnings. These establish scoped computational integrity and traceability, not biological validity or journal acceptance.

## Next work

1. Resolve GSE72833 deposited units, feature identities and the cytology-to-library crosswalk, then specify any held-out evaluation before inspecting its outcomes.
2. Clarify source animals, pooling, batches and phenotype/control details, including N6amt1 age and Rif1 alleles/Cre/crosswalk.
3. Obtain independent scientific review of the conditional models and the publication contribution. Avoid selecting the favorable correlation policy as confirmation.

## Reproduce the presentation

```bash
bash envs/python.sh scripts/plot_research_progress.py
bash envs/python.sh scripts/check_all_snapshots.py
```

The figures read frozen public snapshots. Their stage descriptions are documentary assessments, not inferred quality scores. Original analysis code, contracts and snapshots remain unchanged. The unpublished manuscript, research registries, raw payloads and local proxy configuration stay outside the public release.
