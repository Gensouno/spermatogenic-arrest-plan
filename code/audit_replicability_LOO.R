#!/usr/bin/env Rscript
# Verify all 40 count-cohort one-sample omissions against the prior cache.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
out<-'results/replicability';dir.create(out,recursive=TRUE,showWarnings=FALSE)
cohorts<-c('GSE236388','GSE283904','GSE68283','GSE68284','GSE101511')
old<-fread('results/robustness/leave_one_out_logFC.csv.gz')
old<-old[gse %in% cohorts]
stopifnot(!anyDuplicated(old[,.(gse,variant,gene)]))
records<-list();audit<-list();index<-0L
for(cohort in cohorts) {
  d<-read.delim(gzfile(file.path('data/bulk',paste0(cohort,'_matrix.tsv.gz'))),check.names=FALSE)
  x<-as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double'
  meta<-fread(file.path('data/bulk',paste0(cohort,'_samples.tsv')))
  stopifnot(identical(colnames(x),meta$matrix_sample),all(is.finite(x)),all(x>=0),all(abs(x-round(x))<1e-6))
  variants<-paste0('drop_',meta$matrix_sample)
  stopifnot(setequal(variants,unique(old[gse==cohort]$variant)))
  for(i in seq_len(ncol(x))) {
    group<-factor(meta$group[-i],levels=c('WT','KO'));design<-model.matrix(~group)
    stopifnot(min(table(group))>=2)
    y<-DGEList(counts=x[,-i,drop=FALSE],group=group)
    keep<-filterByExpr(y,design=design);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
    fit<-eBayes(lmFit(voom(y,design,plot=FALSE),design))
    effect<-fit$coefficients[,'groupKO'];variant_name<-variants[i]
    previous<-old[gse==cohort & variant==variant_name]
    stopifnot(setequal(names(effect),previous$gene))
    previous<-previous[match(names(effect),gene)]
    error<-max(abs(effect-previous$logFC))
    stopifnot(is.finite(error),error<1e-8)
    index<-index+1L
    records[[index]]<-data.table(gse=cohort,variant=variant_name,gene=names(effect),logFC=as.numeric(effect))
    audit[[index]]<-data.table(gse=cohort,variant=variant_name,dropped=meta$matrix_sample[i],
       tested_genes=length(effect),minimum_group_size=min(table(group)),max_effect_error=error)
  }
  cat('Recomputed and verified count omissions:',cohort,ncol(x),'variants\n')
}
stopifnot(index==40L)
fwrite(rbindlist(records),file.path(out,'recomputed_count_LOO_effects.csv.gz'),compress='gzip')
fwrite(rbindlist(audit),file.path(out,'LOO_recomputation_audit.csv'))
