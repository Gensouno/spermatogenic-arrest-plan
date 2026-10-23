#!/usr/bin/env Rscript
# Optional reviewed expansion: never overwrite the frozen baseline outputs.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)
  library(AnnotationDbi);library(org.Mm.eg.db)})
out<-'results/data_expansion'
dir.create(file.path(out,'de'),recursive=TRUE,showWarnings=FALSE)
dir.create('data/bulk/expansion',recursive=TRUE,showWarnings=FALSE)
cohorts<-c('GSE259242_P9','GSE117706_P2','GSE117706_P7','GSE235777_P45')
summaries<-list();sensitivities<-list();sensitivity_effects<-list()
fit_counts<-function(x,group) {
  design<-model.matrix(~factor(group,levels=c('Control','KO')))
  d<-DGEList(counts=x)
  use<-filterByExpr(d,design=design)
  d<-calcNormFactors(d[use,,keep.lib.sizes=FALSE])
  v<-voom(d,design,plot=FALSE)
  f<-eBayes(lmFit(v,design))
  list(table=topTable(f,coef=2,number=Inf,sort.by='none'),expression=v$E)
}
for(c in cohorts) {
  raw<-read.delim(gzfile(file.path('data/interim/expansion',paste0(c,'_raw.tsv.gz'))),check.names=FALSE)
  ids<-raw[[1]];x<-as.matrix(raw[-1]);storage.mode(x)<-'double'
  meta<-fread(file.path('data/interim/expansion',paste0(c,'_samples.tsv')))
  stopifnot(identical(colnames(x),meta$sample_id),all(is.finite(x)),all(x>=0),all(abs(x-round(x))<1e-6))
  keys<-sub('\\.[0-9]+$','',ids)
  mapped<-mapIds(org.Mm.eg.db,keys=unique(keys),keytype='ENSEMBL',column='SYMBOL',multiVals='asNA')
  symbols<-unname(mapped[keys]);use<-!is.na(symbols)&nzchar(symbols)
  counts<-rowsum(x[use,,drop=FALSE],symbols[use],reorder=FALSE)
  mass<-colSums(counts)/colSums(x)
  if(min(mass)<.95) stop('Poor mapping coverage: ',c)
  fwrite(data.frame(input_id=ids,symbol=symbols,retained=use),file.path(out,paste0(c,'_mapping.csv')))
  f<-fit_counts(counts,meta$group);tt<-f$table
  fwrite(data.frame(gene=rownames(tt),tt,row.names=NULL),file.path(out,'de',paste0(c,'_de.csv')))
  con<-gzfile(file.path('data/bulk/expansion',paste0(c,'_matrix.tsv.gz')),'wt')
  write.table(data.frame(gene=rownames(counts),counts,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE)
  close(con)
  fwrite(meta,file.path(out,paste0(c,'_sample_evidence.csv')))
  pca<-prcomp(t(f$expression),center=TRUE,scale.=FALSE);v<-pca$sdev^2/sum(pca$sdev^2)
  fwrite(data.frame(gsm=meta$sample_id,group=meta$group,pca$x[,1:2],PC1_variance=v[1],PC2_variance=v[2]),
         file.path(out,paste0(c,'_pca.csv')))
  summaries[[c]]<-data.frame(cohort=c,gse=meta$gse[1],age=meta$age[1],role=meta$role[1],
    control_genotype=paste(unique(meta$genotype[meta$group=='Control']),collapse=';'),
    KO_genotype=paste(unique(meta$genotype[meta$group=='KO']),collapse=';'),Control=sum(meta$group=='Control'),KO=sum(meta$group=='KO'),
    symbol_genes=nrow(counts),tested_genes=nrow(tt),DEG_FDR05_absLFC1=sum(tt$adj.P.Val<.05&abs(tt$logFC)>1),
    min_retained_count_fraction=min(mass),contrast='KO minus explicitly stated Control')
  # Leave-one-out diagnostics, with no automatic exclusions and no CI claims.
  for(j in seq_len(ncol(counts))) {
    refit<-fit_counts(counts[,-j,drop=FALSE],meta$group[-j])$table
    common<-intersect(rownames(tt),rownames(refit))
    sensitivities[[paste(c,j)]]<-data.frame(cohort=c,deleted=meta$sample_id[j],genes=length(common),
      spearman_with_full=cor(tt[common,'logFC'],refit[common,'logFC'],method='spearman'))
    sensitivity_effects[[paste(c,j)]]<-data.frame(cohort=c,deleted=meta$sample_id[j],
      gene=rownames(refit),logFC=refit$logFC)
  }
  print(summaries[[c]])
}
fwrite(rbindlist(summaries),file.path(out,'cohort_summary.csv'))
fwrite(rbindlist(sensitivities),file.path(out,'sample_sensitivity.csv'))
con<-gzfile(file.path(out,'leave_one_out_logFC.csv.gz'),'wt')
write.table(rbindlist(sensitivity_effects),con,sep=',',quote=FALSE,row.names=FALSE)
close(con)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
