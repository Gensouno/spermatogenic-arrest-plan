#!/usr/bin/env Rscript
# Four local gene-count cohorts. Additional transcript/TPM/array cohorts need
# separate importers and are not silently treated as comparable count matrices.
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})
dir.create('data/bulk', recursive=TRUE, showWarnings=FALSE)
dir.create('results/de', recursive=TRUE, showWarnings=FALSE)
dir.create('results/import_qc', recursive=TRUE, showWarnings=FALSE)
sheet <- fread('data/samples.tsv')
files <- c(GSE236388='GSE236388_XAP5_counts.txt.gz',
           GSE283904='GSE283904_Gene_count.txt.gz',
           GSE68283='GSE68283_Raw_counts.txt.gz',
           GSE68284='GSE68284_Raw_counts.txt.gz')
audits <- list()
for (cohort in names(files)) {
  gse <- cohort
  f <- file.path('data/raw', gse, files[[gse]])
  # Base gzfile supports gzip without depending on an optional fread helper.
  tab <- read.delim(gzfile(f), check.names=FALSE, stringsAsFactors=FALSE)
  ids <- as.character(tab[[1]])
  x <- as.matrix(tab[-1]); storage.mode(x) <- 'double'
  stopifnot(!anyNA(x), all(is.finite(x)), all(x >= 0),
            all(abs(x-round(x)) < 1e-6), !anyDuplicated(colnames(x)))
  samples <- sheet[gse == cohort & tissue_class == 'whole_testis' & group %in% c('WT','KO')]
  # Explicit, auditable aliasing between matrix column labels and GSM titles.
  aliases <- samples$title
  if (gse == 'GSE236388') {
    aliases <- ifelse(samples$group == 'WT',
                      paste0('XAP5 WT replicate', sub('.*rep', '', samples$title)),
                      paste0('XAP5 cKO replicate', sub('.*rep', '', samples$title)))
  }
  if (gse == 'GSE283904') {
    aliases <- ifelse(samples$group == 'WT',
                      paste0('N6_CTL_', sub('.*-', '', samples$title)),
                      paste0('N6_KO_', sub('.*-', '', samples$title)))
  }
  samples[, matrix_sample := aliases]
  stopifnot(!anyDuplicated(aliases), setequal(aliases, colnames(x)))
  samples <- samples[match(colnames(x), matrix_sample)]
  stopifnot(identical(samples$matrix_sample, colnames(x)),
            all(table(samples$group) >= 3), !anyDuplicated(samples$sample_id))
  symbols <- ids
  if (all(grepl('^ENSMUSG', ids))) {
    keys <- sub('\\.[0-9]+$', '', ids)
    mapped <- mapIds(org.Mm.eg.db, keys=unique(keys), keytype='ENSEMBL',
                     column='SYMBOL', multiVals='asNA')
    symbols <- unname(mapped[keys])
  }
  keep <- !is.na(symbols) & nzchar(symbols)
  mapping <- data.frame(input_id=ids, symbol=symbols, retained=keep)
  fwrite(mapping, file.path('results/import_qc', paste0(gse, '_gene_mapping.csv')))
  # Multiple original rows assigned to one symbol are summed, never duplicated.
  counts <- rowsum(x[keep,,drop=FALSE], symbols[keep], reorder=FALSE)
  stopifnot(all(colSums(counts) > 0))
  connection <- gzfile(file.path('data/bulk', paste0(gse, '_matrix.tsv.gz')), 'wt')
  write.table(data.frame(gene=rownames(counts), counts, check.names=FALSE), connection,
              sep='\t', quote=FALSE, row.names=FALSE)
  close(connection)
  fwrite(samples, file.path('data/bulk', paste0(gse, '_samples.tsv')), sep='\t')
  group <- factor(samples$group, levels=c('WT','KO'))
  design <- model.matrix(~group)
  d <- DGEList(counts=counts, group=group)
  expressed <- filterByExpr(d, design=design)
  d <- calcNormFactors(d[expressed,,keep.lib.sizes=FALSE], method='TMM')
  v <- voom(d, design=design, plot=FALSE)
  fit <- eBayes(lmFit(v, design))
  result <- topTable(fit, coef='groupKO', number=Inf, sort.by='none')
  result <- data.frame(gene=rownames(result), result, row.names=NULL)
  fwrite(result, file.path('results/de', paste0(gse, '_de.csv')))
  pca <- prcomp(t(v$E), center=TRUE, scale.=FALSE)
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  scores <- data.frame(sample=colnames(counts), group=group, pca$x[,1:2],
                       PC1_variance=variance[1], PC2_variance=variance[2])
  fwrite(scores, file.path('results/import_qc', paste0(gse, '_pca.csv')))
  audits[[gse]] <- data.frame(gse=gse, n_samples=ncol(counts),
      WT=sum(group=='WT'), KO=sum(group=='KO'), input_rows=nrow(x),
      unmapped_or_ambiguous=sum(!keep), symbol_genes=nrow(counts), tested_genes=nrow(result),
      DEG_FDR05_absLFC1=sum(result$adj.P.Val < .05 & abs(result$logFC)>1),
      comparison='KO minus WT', normalization='TMM + voom',
      grouping_source='data/samples.tsv; see sample_metadata stage official_sample_evidence.csv for GEO source evidence')
  print(audits[[gse]])
}
fwrite(rbindlist(audits), 'results/import_qc/count_cohorts.csv')
capture.output(sessionInfo(), file='results/import_qc/R_sessionInfo.txt')
