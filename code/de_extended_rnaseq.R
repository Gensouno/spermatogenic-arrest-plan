#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
  library(AnnotationDbi); library(org.Mm.eg.db)
})
audits <- list()
for (gse in c('GSE101511','GSE294209')) {
  tab <- read.delim(gzfile(file.path('data/interim/extended',paste0(gse,'_raw.tsv.gz'))),check.names=FALSE)
  ids <- as.character(tab[[1]])
  x <- as.matrix(tab[-1]); storage.mode(x) <- 'double'
  if (gse == 'GSE101511') {
    summary_rows <- startsWith(ids,'__')
    fwrite(data.frame(category=ids[summary_rows],x[summary_rows,,drop=FALSE],check.names=FALSE),
           'results/import_qc/GSE101511_HTSeq_summary.csv')
    x <- x[!summary_rows,,drop=FALSE]
    ids <- ids[!summary_rows]
  }
  meta <- fread(file.path('data/bulk',paste0(gse,'_samples.tsv')))
  stopifnot(identical(colnames(x),meta$matrix_sample), all(is.finite(x)), all(x>=0))
  symbols <- ids
  if (gse == 'GSE101511') {
    keys <- sub('\\.[0-9]+$','',ids)
    legacy <- fread('data/ref/annotation/ensembl79_transcript_mapping.tsv')
    stopifnot(!anyDuplicated(legacy$transcript))
    symbols <- legacy$symbol[match(keys,legacy$transcript)]
  }
  keep <- !is.na(symbols) & nzchar(symbols)
  fwrite(data.frame(input_id=ids,symbol=symbols,retained=keep),
         file.path('results/import_qc',paste0(gse,'_gene_mapping.csv')))
  expr <- rowsum(x[keep,,drop=FALSE], symbols[keep],reorder=FALSE)
  con <- gzfile(file.path('data/bulk',paste0(gse,'_matrix.tsv.gz')),'wt')
  write.table(data.frame(gene=rownames(expr),expr,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE)
  close(con)
  group <- factor(meta$group,levels=c('WT','KO'))
  design <- model.matrix(~group)
  if (gse == 'GSE101511') {
    d <- DGEList(counts=expr,group=group)
    expressed <- filterByExpr(d,design=design)
    d <- calcNormFactors(d[expressed,,keep.lib.sizes=FALSE])
    v <- voom(d,design,plot=FALSE)
    e <- v$E
    fit <- eBayes(lmFit(v,design))
    method <- 'Ensembl79 transcript features summed by symbol; TMM+voom'
    unit <- 'counts'
  } else {
    # Author-supplied TPM is a compositional abundance, never raw counts.
    expressed <- rowSums(expr>=1)>=3
    e <- log2(expr[expressed,,drop=FALSE]+0.5)
    fit <- eBayes(lmFit(e,design),trend=TRUE,robust=TRUE)
    method <- 'log2(TPM+0.5), limma-trend; exploratory sensitivity'
    unit <- 'TPM'
  }
  tt <- topTable(fit,coef='groupKO',number=Inf,sort.by='none')
  fwrite(data.frame(gene=rownames(tt),tt,row.names=NULL),file.path('results/de',paste0(gse,'_de.csv')))
  pca <- prcomp(t(e),center=TRUE,scale.=FALSE)
  frac <- pca$sdev^2/sum(pca$sdev^2)
  fwrite(data.frame(sample=colnames(expr),group=group,pca$x[,1:2],PC1_variance=frac[1],PC2_variance=frac[2]),
         file.path('results/import_qc',paste0(gse,'_pca.csv')))
  retained_fraction <- colSums(x[keep,,drop=FALSE])/colSums(x)
  fwrite(data.frame(sample=colnames(x),original_total=colSums(x),retained_total=colSums(expr),retained_fraction),
         file.path('results/import_qc',paste0(gse,'_mapping_mass.csv')))
  audits[[gse]] <- data.frame(gse=gse,unit=unit,n_samples=ncol(expr),WT=sum(group=='WT'),KO=sum(group=='KO'),
      input_rows=nrow(x),unmapped_or_ambiguous=sum(!keep),symbol_genes=nrow(expr),tested_genes=nrow(tt),
      min_retained_abundance_fraction=min(retained_fraction),
      DEG_FDR05_absLFC1=sum(tt$adj.P.Val<.05 & abs(tt$logFC)>1),method=method)
  print(audits[[gse]])
}
fwrite(rbindlist(audits),'results/import_qc/extended_rnaseq_cohorts.csv')
capture.output(sessionInfo(),file='results/import_qc/extended_R_sessionInfo.txt')
