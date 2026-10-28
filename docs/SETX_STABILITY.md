# Setx count-model omission stability and matched-universe diagnostics

This is a post-hoc follow-up on the same GSE224875 libraries. It completes the count-model omission analysis left open in the preceding [source and length report](SETX_SOURCE_AND_LENGTH.md). It does not supply new biological replication, resolve source-code/animal uncertainties, or change the original four-unit analysis's zero adjusted directional discoveries. The 105 descriptive candidates and 36 downward priorities remain fixed.

## Prespecified scope and frozen inputs

The new analysis contract is `data/setx_stability_contract.json`, specified on 2026-10-28 after the preceding TPM/count outcomes had been inspected. `data/setx_stability_inputs.json` binds the original RSEM TAR, archived sample/member identities, gene mapping, candidate lists and preceding fit outputs. Each checksum was checked against earlier frozen identities before being recorded. Earlier sources and snapshots were not rewritten.

Original `expected_count`, `effective_length` and `TPM` values were reconstructed from all 40 genome members. Estimated counts were not rounded. The frozen native edgeR length-offset fitting function was reused with expression filtering, positive-length checks, TMM normalization, robust quasi-likelihood fitting and a four-column categorical age/genotype interaction contrast. Its full and ambiguous-pair-exclusion fits reproduce the preceding values to within `1e-9` maximum absolute error.

## Count omission suite

The suite contains **42 fits**: all 40 libraries, exclusion of both unresolved P16 records, and every one of the 40 individual library omissions. Filtering and normalization are repeated after each omission. Every variant adjusts its complete tested original-gene family with BY and BH; candidate annotations are applied afterwards. A missing, filtered-out or non-uniquely mapped candidate fails retention. Omissions are sensitivity checks, not independent animals or confidence intervals.

For the fixed 36 priorities, all 36 remain measured across the full suite:

| Descriptive property across all 42 count fits | Priorities / 36 |
|---|---:|
| Negative P18-minus-P10 change in the mutant-minus-control contrast | 35 |
| Change at most `-log2(1.5)` in every fit | 33 |
| Conditional BY value below 0.05 in every fit's complete gene family | 13 |
| Magnitude criterion retained in both this suite and the prior TPM suite | 31 |

The last two rows are descriptive intersections. They are not new FDR-controlled candidate sets, tests that the true interaction exceeds a magnitude threshold, or evidence of a common causal mechanism.

**Adcy1** crosses zero across omissions. **Casz1** and **Fsip2** remain negative but do not meet the magnitude threshold in every variant. These findings are retained without changing labels or excluding inconvenient libraries. The complete 105-candidate table records missingness and the variant giving the least negative estimated change.

## Matched-universe diagnostics

The preceding full TPM and count models shared **16,853 original gene IDs**. This intersection is fixed before the new diagnostic fits. All three diagnostics use those same genes, the same 40 libraries and the same categorical design, with no additional gene removal.

1. Refit log2(TPM + 0.5) with trend/robust limma on the common universe.
2. Refit the estimated-count model with centered gene/sample length offsets.
3. Refit the same count matrix with those offsets set to zero, as an intentionally incomplete counterfactual.

The two count diagnostics share exact library sizes and TMM factors. Their fitting policy, coefficient prior count and contrast are identical, while dispersions are reestimated under each offset. The numerical audit finds no difference in their library-normalization component. A synthetic eightfold count-and-length shift verifies the expected effect change when the offset is toggled; constant per-gene lengths produce numerically equivalent fits.

The zero-offset diagnostic is **not a recommended alternative workflow**. The [official tximport documentation](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html) cautions against passing uncorrected original gene estimates downstream without appropriate length correction. Here the deliberately changed specification probes model dependence; it does not certify a preferred biological answer.

| Fixed-universe diagnostic | Tested genes | Conditional BY<0.05, all genes | Fixed priorities with conditional BY<0.05 |
|---|---:|---:|---:|
| TPM, limma | 16,853 | 8,419 | 32 / 36 |
| Count, length offset | 16,853 | 5,857 | 17 / 36 |
| Count, zero-offset counterfactual | 16,853 | 6,745 | 20 / 36 |

Thus the candidate-support difference persists with an identical hypothesis universe; it cannot be explained solely by different gene coverage or multiplicity-family sizes. Turning off the length prior does not make the count model equivalent to TPM-limma. Offset presence affects both fitted means and subsequent dispersion estimates, and TPM/count scale, normalization and mean-variance modeling still differ. These diagnostics do **not** assign a causal percentage of the original difference to length correction.

Threshold memberships are not nested: Drc7 meets the conditional BY threshold only with the length offset among the two count diagnostics, while Casz1, Cyp2u1, Fsip2 and Spef2 meet it only in the zero-offset counterfactual. Net rejection-count changes must not be described as an exact list of genes simply gained or lost. No candidate-only adjustment is used.

![Fixed Setx priorities under count omissions and matched models. Left: the full count estimate, its range over 42 count fits, and the prior full TPM estimate; the range is not a confidence interval. Right: conditional BY threshold status in three newly fitted models on the same gene universe. The zero-offset column is a diagnostic counterfactual.](../setx_stability_snapshot/stability_and_matched_models.png)

The left panel uses each original complete-gene-family fit; the right panel uses the matched-universe refits. Their denominators and roles are intentionally distinct. All 36 priorities are displayed, including model-sensitive and nonsignificant entries.

## Execution and verification

```bash
# Reconstruct the original profiles; run the full omission suite and diagnostics.
bash envs/python.sh scripts/run_setx_stability.py

# Check complete test families and fixed normalization; reproduce in a fresh directory.
bash envs/python.sh scripts/check_setx_stability.py --reproduce
```

Independent checks recompute all 42 omission-family and three matched-family BY/BH adjustments, verify original-gene inventories, sample counts, normalization equality and agreement with the preceding frozen fits. Fresh-directory reproduction rereads the original processed TAR, refits all models, and compares **10 tables, 2 JSON records and 1 figure**, plus the three complete reconstructed matrices. It reuses installed dependencies and frozen annotation/prior comparisons; there is no fresh network test, raw-read alignment, clean installation or biological validation.

The compact snapshot retains candidate-level omissions, model summaries, complete matched-universe statistics, normalization audits, the figure and reproduction record. The larger complete-gene omission table remains in ignored `results/setx_stability/` and is included in the numerical reproduction comparison.

## Interpretation and next research priorities

A substantial subset retains a descriptive magnitude under both expression approaches, but statistical support is model-dependent. Such agreement is useful for prioritizing follow-up, not a substitute for source animal/pooling/technical verification or independent evidence. Code `7844`, barcode nesting, source independence and quantification uncertainty remain unresolved.

This closes the immediate same-source omission question. Further work should prioritize independent age-compatible data and scientific review rather than expand the same-source variant count without a new methodological reason. No new shared cell-intrinsic mechanism, calibrated fraction estimate or clinical relevance is established.
