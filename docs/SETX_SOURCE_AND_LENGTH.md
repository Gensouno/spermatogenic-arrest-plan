# Setx source audit and length-aware count sensitivity

This follow-up audits the identities behind GSE224875 and checks whether the previous exploratory age-by-genotype observations depend on TPM-only modeling. It uses the same source libraries, not an independent replication. The original four-genetic-unit directional family still has zero adjusted discoveries. Earlier snapshots and their bound sources are preserved.

## Retrieval and sample identity

The targeted retrieval contract is `data/setx_source_audit_contract.json`. Sources were accessed on 2026-10-27 through the configured local proxy. Current GEO SOFT was compared with the earlier frozen record; ENA return/search fields were inspected, the project run count was queried first, and all 40 records were retrieved with an explicit limit of 100. Exact URLs, hashes and access failures are retained in the retrieval audit; no direct-network fallback was used.

All 40 frozen GSM records reconcile to 40 current GEO records, 40 BioSamples, 40 SRA experiments and 40 ENA runs. Every run is reported as paired-end RNA-Seq on an Illumina HiSeq 2500; a shared instrument model is not proof of a shared physical instrument or batch. Titles, genotype labels and supplementary-file identities agree with the frozen source. Distinct database accessions do not establish distinct animals.

The two P16 libraries with filename code `7844` remain unresolved:

| Source library | Source label | BioSample | Experiment | Run |
|---|---|---|---|---|
| GSM7034763 | WT P16 | SAMN33213895 | SRX19315358 | SRR23374748 |
| GSM7034783 | SetxKO P16 | SAMN33213875 | SRX19315378 | SRR23374728 |

Their separate sample, experiment and run XML records preserve the original labels but do not supply an animal/pooling manifest or clarify the filename code. No relabeling is justified. Barcode nesting identified in the preceding branch remains unresolved; the audit does not establish the absence or presence of a technical artifact.

## Allele, phenotype and publication versions

The GEO summary identifies **spcar3**. The matching publication is [Fujiwara et al., Reproduction 2023, DOI 10.1530/REP-23-0166](https://pubmed.ncbi.nlm.nih.gov/37801077/). Its abstract reports meiotic arrest and a lack of spermatids/mature sperm in the mutant. The [version-1 preprint](https://www.biorxiv.org/content/10.1101/2023.04.12.536672v1.full) explicitly links GSE224875 and describes an ENU-induced Setx allele and four biological replicates for P18 RNA-seq. The P18 statement does not establish animal independence at every age or match histology to each reanalyzed library. “KO” is retained as the source's contrast label, not as proof of an engineered knockout.

Publisher HTML/PDF downloads failed; the final abstract was checked in PubMed, while detailed methods were checked in the explicitly versioned preprint HTML/JATS. These are not interchangeable publication versions. The bioRxiv summary API omitted an author and returned `published=NA`, whereas versioned HTML/JATS and Crossref list eight authors and PubMed identifies the published article. Document-level authorship was used, and the inconsistent API summary was retained as an audit limitation.

## Length-aware sensitivity design

The separate contract is `data/setx_length_contract.json`. Original RSEM `expected_count` and `effective_length` were reread from all 40 genome members with frozen member hashes. Estimated counts were never rounded into raw counts. Gene-level effective lengths were reused as supplied; no transcript-level quantification, bootstrap/inferential replicates or quantification-uncertainty correction was reconstructed.

The installed edgeR 4.10.6 provides `DGEListFromTximport`, as documented in the [official tximport workflow](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html). A list of the gene-level estimate and length matrices with `countsFromAbundance="no"` creates centered log-length offset priors. Native offset values were checked independently. Expression filtering, positive-length checks, TMM library normalization and robust edgeR v4 quasi-likelihood fitting were applied with ten categorical age/genotype cell means. The same four interaction contrasts as the preceding TPM branch were tested jointly. BY adjustment covers every tested original gene ID in each variant; BH is a sensitivity. No candidate-only adjustment was performed.

Two variants were fitted: all 40 libraries and the 38-library sensitivity excluding both unresolved P16 records. No expressed genes were removed for nonpositive length in either observed-data fit. The default coefficient prior count was explicitly fixed at 0.125; these count-model coefficients and the previous differences of mean log-TPM are different estimands. This branch does not run a count-based single-library omission suite and cannot inherit the prior branch's all-omission stability claim.

## Results and model dependence

| Model | Tested original gene IDs | Conditional BY<0.05 | Conditional BH<0.05 |
|---|---:|---:|---:|
| Length offset, 40 libraries | 18,207 | 5,891 | 9,098 |
| Length offset, exclude P16 pair | 18,470 | 6,599 | 9,515 |

The 40-library count and preceding TPM fits share 16,853 tested gene IDs. Within that matched universe, conditional BY thresholds are met by 5,843 count-model tests, 8,345 TPM tests, and 4,667 under both. Neither overlapping rejection counts nor the filtering differences define an additional calibrated discovery family.

All 36 fixed downward priorities are tested and uniquely mapped in both length-offset fits. In each fit, all 36 have a negative P18-minus-P10 change in genotype contrast and 34 meet the descriptive magnitude threshold of at most `-log2(1.5)`. In the full count fit, 35/36 signs agree with the full TPM fit; **Adcy1** changes from a small positive to a small negative contrast, and neither model rejects its omnibus test at the BY threshold. Conditional BY rejections among the fixed priorities number 17/36 in the full count fit and 24/36 after pair exclusion. These are annotations of the full gene families, not subgroup FDR guarantees.

Thus many effect directions persist under a length-aware model, while inferential strength and some effects depend on the expression unit, filtering, normalization and model. These results do not establish intrinsic regulation, validated cell fractions or a common causal mechanism. The same unresolved biological and technical design applies to all variants.

## Reproduction and checks

```bash
# Offline audit of frozen GEO/ENA/publication responses.
bash envs/python.sh scripts/audit_setx_sources.py

# Rebuild the RSEM estimate/length matrices and fit the two QL models.
bash envs/python.sh scripts/run_setx_length.py

# Independently recompute complete-family adjustments and reproduce in a new directory.
bash envs/python.sh scripts/check_setx_followup.py --reproduce
```

Fresh-directory reproduction matches **8 tables and 4 JSON records**, plus the complete reconstructed count/length matrices. It rereads the frozen source TAR and metadata, refits both models, and reuses installed dependencies, prior annotation and the frozen TPM comparison. It does not repeat network downloads, raw-read alignment, a clean installation or biological experiments.

Five new regression tests cover identity/count mismatches, distinct-accession versus animal-independence semantics, ambiguous source relations, fractional-count preservation, nonpositive-length exclusion, length offsets and a known interaction direction. Existing scientific snapshots remain unmodified. The next priorities are source clarification of animal/pooling/technical design and independent validation; successful arithmetic and source reconciliation do not resolve those questions.

## Software workflow attribution

The database retrieval used the installed `database-lookup` skill. Its procedural source is Kassis, T.; Agarwal, V.; He, Y.; Patel, D.; Brueckner, A. M. (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*, [arXiv:2609.00065](https://arxiv.org/abs/2609.00065). The current arXiv metadata was checked; this attribution does not validate the statistical or biological conclusions.
