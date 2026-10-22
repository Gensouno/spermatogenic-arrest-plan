#!/usr/bin/env Rscript
# HISTORICAL IMPLEMENTATION, not the default pipeline. Gene reconstruction
# R2 does not validate cell fractions; see docs/METHODS.md and docs/RESULTS.md.
# =============================================================================
# 05_deconvolution.R -- estimate cell-type proportions P for every bulk sample
#
# Model:  X[g,s]  ~  sum_k S[g,k] * P[k,s]
#   X = bulk expression (genes x samples), S = cell-type signature (genes x K)
#   P = proportions (K x samples), the unknown we solve for.
#
# We solve per sample by non-negative least squares, then rescale so the
# proportions sum to 1. NNLS is used instead of MuSiC/BisqueRNA because those
# two are GitHub-only (verified NOT on CRAN or Bioconductor 3.23) and the
# restricted network access makes GitHub installs unreliable; nnls + quadprog are on
# CRAN and install cleanly into the project-local .Rlib.
#
# CRITICAL assumption being tested by this stage: deconvolution of BULK WHOLE
# TESTIS is only meaningful when the sample really is whole testis. A sorted
# pachytene-spermatocyte population has its composition fixed by the sort, not
# measured -- for those samples P is not an estimate of biology and the
# residual step is meaningless. The sample sheet must carry a `tissue_class`
# column and this script refuses to deconvolve anything but whole testis.
#
# Usage:
#   Rscript code/05_deconvolution.R
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

.libPaths(c(".Rlib", .libPaths()))

sig_file    <- "data/ref/signature_GSE112393_cpm.csv.gz"
sample_sheet<- "data/samples.tsv"
bulk_dir    <- "data/bulk"
out_dir     <- "results/deconvolution"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(sig_file), file.exists(sample_sheet))

if (!requireNamespace("nnls", quietly = TRUE))
  stop("nnls missing -- run envs/setup.sh (installs into .Rlib)")

# -----------------------------------------------------------------------------
# Signature matrix S: rows = cell types, cols = genes
# The builder writes a CSV with the cell type in the first column.
# -----------------------------------------------------------------------------
Sdt <- fread(sig_file)
ct <- Sdt[[1]]
Sdt[, (1) := NULL]
S <- as.matrix(Sdt)
rownames(S) <- ct
cat("signature:", nrow(S), "cell types x", ncol(S), "genes\n")
cat("cell types:", paste(rownames(S), collapse = ", "), "\n")

# Deconvolution wants the signature on a linear (not log) scale. The builder
# already writes mean counts per cell, so we only need to drop all-zero genes.
S <- S[rowSums(S) > 0, , drop = FALSE]
cat("after dropping all-zero genes:", nrow(S), "x", ncol(S), "\n")

samples <- fread(sample_sheet)
if (!"tissue_class" %in% names(samples))
  stop("data/samples.tsv must have a tissue_class column (whole_testis / sorted / scrna)")

targets <- samples[tissue_class == "whole_testis"]
if (nrow(targets) == 0)
  stop("no whole_testis samples to deconvolve -- check samples.tsv")
skipped <- samples[tissue_class != "whole_testis"]
cat(sprintf("deconvolving %d whole-testis samples; skipping %d non-whole-testis\n",
            nrow(targets), nrow(skipped)))
if (nrow(skipped)) print(skipped[, .(gse, tissue_class, n = .N), by = .(gse, tissue_class)])

# -----------------------------------------------------------------------------
# NNLS per sample, with a held-out-gene cross-validation for the fit quality.
#
# Why cross-validate: a proportion vector can be produced for ANY sample no
# matter how badly the reference matches. Without an out-of-sample error we
# cannot tell "composition really is X" from "the reference cannot represent
# this tissue". We fit on 80% of genes and score on the held-out 20%, which
# gives an honest R^2 per sample. Samples with low held-out R^2 are flagged and
# excluded from the composition null model (their P is not trustworthy).
# -----------------------------------------------------------------------------
set.seed(42)

fit_nnls <- function(x, S) {
  # x: gene vector (named), S: K x G matrix, both on the same gene order
  A <- t(S)                       # G x K, the design
  sol <- nnls::nnls(A, x)
  p <- sol$x
  if (sum(p) > 0) p <- p / sum(p)
  names(p) <- rownames(S)
  list(p = p, rss = sum((A %*% sol$x - x)^2), coef = sol$x)
}

deconvolve_one <- function(mat, S, genes) {
  # mat: genes x samples
  G <- length(genes); idx <- seq_len(G)
  tr <- sample(idx, floor(0.8 * G))
  te <- setdiff(idx, tr)
  out <- list()
  for (s in colnames(mat)) {
    x <- as.numeric(mat[, s]); names(x) <- genes
    # full-data fit -> the proportions we report
    f <- fit_nnls(x, S)
    # train/test fit -> honest R^2
    Sr <- S[, tr, drop = FALSE]
    ft <- fit_nnls(x[tr], Sr)
    pred <- as.numeric(t(S[, te, drop = FALSE]) %*% ft$coef)
    obs  <- x[te]
    r2 <- 1 - sum((obs - pred)^2) / max(sum((obs - mean(obs))^2), 1e-12)
    out[[s]] <- list(p = f$p, r2_heldout = r2,
                     n_genes = G, n_cells_sig = NA)
  }
  out
}

all_P <- list(); qc <- list()
for (gse in unique(targets$gse)) {
  f <- file.path(bulk_dir, paste0(gse, "_matrix.tsv.gz"))
  if (!file.exists(f)) { cat("MISSING", f, "-- skipped\n"); next }
  M <- as.matrix(fread(f))
  genes <- rownames(M)
  common <- intersect(genes, colnames(S))
  if (length(common) < 2000)
    warning(sprintf("%s: only %d genes shared with signature", gse, length(common)))
  M <- M[common, , drop = FALSE]
  # Signature on the same gene order and same linear scale.
  Sg <- S[, common, drop = FALSE]
  cat(sprintf("%s: %d genes x %d samples\n", gse, nrow(M), ncol(M)))
  res <- deconvolve_one(M, Sg, common)
  for (s in names(res)) {
    all_P[[s]] <- res[[s]]$p
    qc[[s]] <- data.table(sample = s, gse = gse,
                          r2_heldout = res[[s]]$r2_heldout,
                          n_genes = res[[s]]$n_genes)
  }
}

P <- do.call(cbind, all_P)
Q <- rbindlist(qc)
P <- P[, Q$sample, drop = FALSE]     # align columns to the QC table

saveRDS(list(P = P, S = S, qc = Q), file.path(out_dir, "deconvolution.rds"))

fwrite(as.data.table(P, keep.rownames = "celltype"),
       file.path(out_dir, "proportions.csv"))
fwrite(Q, file.path(out_dir, "deconvolution_qc.csv"))

cat("\n=== proportion summary by cell type ===\n")
print(round(rowMeans(P), 4))
cat("\n=== held-out R^2 distribution ===\n")
print(summary(Q$r2_heldout))
bad <- Q[r2_heldout < 0.5]
if (nrow(bad)) {
  cat("\nWARNING: samples with held-out R^2 < 0.5 (P not trustworthy):\n")
  print(bad)
  fwrite(bad, file.path(out_dir, "flagged_samples.csv"))
}
cat("\nwrote:", file.path(out_dir, "proportions.csv"),
    "/ deconvolution_qc.csv\n")

# -----------------------------------------------------------------------------
# Diagnostic figure: stacked composition per sample, KO vs WT side by side.
# This is the raw material for Act II -- if KO and WT bars look different, there
# IS a composition shift, and the question becomes how much of the expression
# similarity it explains rather than whether it exists at all.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2) })
Pd <- as.data.table(P, keep.rownames = "celltype")
long <- melt(Pd, id.vars = "celltype", variable.name = "sample", value.name = "frac")
long <- merge(long, samples[, .(sample, gse, group, tissue_class)],
              by = "sample", all.x = TRUE)
long <- long[!is.na(group)]
p <- ggplot(long, aes(sample, frac, fill = celltype)) +
  geom_col(width = 0.85) +
  facet_grid(. ~ gse, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "estimated fraction", fill = NULL,
       title = "Deconvolved cell-type composition (whole testis only)",
       subtitle = paste0("GSE112393 signature; ", nrow(Q), " samples")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "bottom")
ggsave("figures/05_composition_stacked.png", p, width = 12, height = 5, dpi = 150)
cat("wrote: figures/05_composition_stacked.png\n")
