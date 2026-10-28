#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma);library(AnnotationDbi);library(org.Mm.eg.db)})
out<-commandArgs(trailingOnly=TRUE)[1];meta<-fread(file.path(out,'samples.csv'))
priorities<-fread('replicability_snapshot/stable_descriptive_priorities.csv')$gene
allstats<-list();effects<-list();summaries<-list();pcas<-list();annotations<-list();libraries<-list()
for(cohort_name in unique(meta$cohort)) {
  mm<-meta[cohort==cohort_name];d<-read.delim(gzfile(file.path(out,paste0(cohort_name,'.tsv.gz'))),check.names=FALSE)
  x<-as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double'
  stopifnot(identical(colnames(x),mm$matrix_sample),ncol(x)==6L,!anyDuplicated(rownames(x)),all(x>=0),all(is.finite(x)),all(x==floor(x)))
  symbol<-if(cohort_name%in%c('Hormad_dsfb_P19','Morc1_P14p5'))unname(mapIds(org.Mm.eg.db,keys=rownames(x),keytype='ENSEMBL',column='SYMBOL',multiVals='asNA')) else rownames(x)
  unique_symbol<-!is.na(symbol)&nzchar(symbol)&!duplicated(symbol)&!duplicated(symbol,fromLast=TRUE)
  annotations[[cohort_name]]<-data.table(cohort=cohort_name,gene=rownames(x),symbol=symbol,globally_unique_symbol=unique_symbol)
  for(i in 0:6) {
    use<-if(i==0)seq_len(6) else setdiff(seq_len(6),i)
    variant<-if(i==0)'full' else paste0('omit_',mm$matrix_sample[i])
    group<-factor(mm$group[use],levels=c('Control','Perturbation'));design<-model.matrix(~group)
    stopifnot(qr(design)$rank==2,min(table(group))>=2)
    y<-DGEList(x[,use,drop=FALSE]);keep<-filterByExpr(y,design=design)
    y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,design,plot=FALSE)
    fit<-eBayes(lmFit(v,design));tt<-topTable(fit,coef=2,number=Inf,sort.by='none')
    tt$BY<-p.adjust(tt$P.Value,'BY');tt$BH<-p.adjust(tt$P.Value,'BH')
    tt$gene<-rownames(tt);tt$symbol<-symbol[match(tt$gene,rownames(x))]
    tt$globally_unique_symbol<-unique_symbol[match(tt$gene,rownames(x))]
    candidate<-tt[tt$globally_unique_symbol&tt$symbol%in%priorities,]
    effects[[paste(cohort_name,variant)]]<-data.table(cohort=cohort_name,variant=variant,candidate)
    summaries[[paste(cohort_name,variant)]]<-data.table(cohort=cohort_name,variant=variant,libraries=length(use),tested_genes=nrow(tt),
      conditional_BY05=sum(tt$BY<.05),conditional_BH05=sum(tt$BH<.05))
    libraries[[paste(cohort_name,variant)]]<-data.table(cohort=cohort_name,variant=variant,matrix_sample=colnames(y),y$samples)
    if(i==0) {
      allstats[[cohort_name]]<-data.table(cohort=cohort_name,tt)
      pc<-prcomp(t(v$E));fraction<-pc$sdev^2/sum(pc$sdev^2)
      pcas[[cohort_name]]<-data.table(cohort=cohort_name,matrix_sample=mm$matrix_sample,group=group,PC1=pc$x[,1],PC2=pc$x[,2],PC1_variance=fraction[1],PC2_variance=fraction[2])
    }
  }
}
fwrite(rbindlist(allstats),file.path(out,'full_statistics.csv.gz'),compress='gzip')
fwrite(rbindlist(effects),file.path(out,'priority_effects.csv.gz'),compress='gzip')
fwrite(rbindlist(summaries),file.path(out,'fit_summary.csv'))
fwrite(rbindlist(pcas),file.path(out,'PCA.csv'))
fwrite(rbindlist(annotations),file.path(out,'gene_annotation.csv.gz'),compress='gzip')
fwrite(rbindlist(libraries),file.path(out,'normalization.csv'))
fwrite(data.table(component=c('R','edgeR','limma','org.Mm.eg.db'),version=c(as.character(getRversion()),vapply(c('edgeR','limma','org.Mm.eg.db'),function(p)as.character(packageVersion(p)),character(1)))),file.path(out,'runtime_R_versions.csv'))
