#!/usr/bin/env Rscript
# Exploratory conditional library-level fits; animal independence not established.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma);library(AnnotationDbi);library(org.Mm.eg.db)})
out<-'results/new_models';prepared<-'data/interim/new_models';meta<-fread(file.path(out,'sample_evidence.csv'))
read_matrix<-function(name){d<-read.delim(gzfile(file.path(prepared,name)),check.names=FALSE);x<-as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double';x}
eef<-read_matrix('Eef1g_counts.tsv.gz');setx<-read_matrix('Setx_TPM.tsv.gz')
annotation<-fread(file.path(out,'Eef1g_source_annotation.csv'))
mapped<-mapIds(org.Mm.eg.db,keys=rownames(setx),keytype='ENSEMBL',column='SYMBOL',multiVals='asNA')
valid<-!is.na(mapped)&nzchar(mapped)
unique_symbol<-valid & !duplicated(mapped) & !duplicated(mapped,fromLast=TRUE)
estimated<-read_matrix('Setx_expected_count.tsv.gz')
stopifnot(identical(rownames(estimated),rownames(setx)),identical(colnames(estimated),colnames(setx)))
fwrite(data.table(gsm=colnames(setx),TPM_sum=colSums(setx),
    named_symbol_TPM_fraction=colSums(setx[valid,,drop=FALSE])/colSums(setx),
    globally_unique_symbol_TPM_fraction=colSums(setx[unique_symbol,,drop=FALSE])/colSums(setx),
    named_symbol_estimated_count_fraction=colSums(estimated[valid,,drop=FALSE])/colSums(estimated),
    globally_unique_symbol_estimated_count_fraction=colSums(estimated[unique_symbol,,drop=FALSE])/colSums(estimated)),
    file.path(out,'source_mapping_mass.csv'))
annotation<-rbindlist(list(data.table(cohort_source='Eef1g',gene=annotation$gene_id,symbol=annotation$symbol),
    data.table(cohort_source='Setx',gene=rownames(setx),symbol=unname(mapped[rownames(setx)]))))
fwrite(annotation,file.path(out,'gene_annotation.csv'))
fit_source<-function(x,mm,kind,pseudo=.5){
  group<-factor(mm$group,levels=c('Control','KO'));design<-model.matrix(~group)
  stopifnot(qr(design)$rank==2,min(table(group))>=2,all(is.finite(x)),all(x>=0))
  if(kind=='counts') {
    y<-DGEList(counts=x);keep<-filterByExpr(y,design=design);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
    v<-voom(y,design,plot=FALSE);f<-eBayes(lmFit(v,design));expression<-v$E
  } else {
    keep<-rowSums(x>=1)>=min(table(group));expression<-log2(x[keep,,drop=FALSE]+pseudo)
    f<-eBayes(lmFit(expression,design),trend=TRUE,robust=TRUE)
  }
  list(table=topTable(f,coef='groupKO',number=Inf,sort.by='none'),expression=expression)
}
tables<-list();loo<-list();summary<-list();pcas<-list();pseudo_rows<-list();index<-0L
for(cohort_name in unique(meta$cohort)) {
  mm<-meta[meta$cohort==cohort_name];kind<-if(startsWith(cohort_name,'Setx'))'TPM' else 'counts'
  stopifnot(uniqueN(mm$cohort)==1L,nrow(mm)==if(kind=='counts')6L else 8L)
  x<-if(kind=='TPM')setx[,mm$matrix_sample,drop=FALSE] else eef[,mm$matrix_sample,drop=FALSE]
  f<-fit_source(x,mm,kind);tt<-f$table
  tables[[cohort_name]]<-data.table(cohort=cohort_name,gene=rownames(tt),tt)
  pc<-prcomp(t(f$expression));v<-pc$sdev^2/sum(pc$sdev^2)
  pcas[[cohort_name]]<-data.table(cohort=cohort_name,gsm=mm$gsm,group=mm$group,PC1=pc$x[,1],PC2=pc$x[,2],PC1_variance=v[1],PC2_variance=v[2])
  summary[[cohort_name]]<-data.table(cohort=cohort_name,unit=kind,samples=ncol(x),Control=sum(mm$group=='Control'),KO=sum(mm$group=='KO'),
    tested_features=nrow(tt),conditional_model_BH05_absLFC1=sum(tt$adj.P.Val<.05&abs(tt$logFC)>1),
    interpretation='Conditional on independent source libraries and model assumptions; no confirmatory biological claim')
  for(i in seq_len(ncol(x))){
    omitted<-fit_source(x[,-i,drop=FALSE],mm[-i],kind)$table;index<-index+1L
    loo[[index]]<-data.table(cohort=cohort_name,deleted=mm$gsm[i],gene=rownames(omitted),logFC=omitted$logFC)
  }
  if(kind=='TPM') for(pseudo in c(.1,1)) {
    alternate<-fit_source(x,mm,kind,pseudo)$table
    pseudo_rows[[paste(cohort_name,pseudo)]]<-data.table(cohort=cohort_name,pseudocount=pseudo,gene=rownames(alternate),logFC=alternate$logFC)
  }
}
fwrite(rbindlist(tables),file.path(out,'conditional_model_statistics.csv.gz'),compress='gzip')
fwrite(rbindlist(loo),file.path(out,'leave_one_out_effects.csv.gz'),compress='gzip')
fwrite(rbindlist(pseudo_rows),file.path(out,'TPM_pseudocount_effects.csv.gz'),compress='gzip')
fwrite(rbindlist(summary),file.path(out,'cohort_summary.csv'));fwrite(rbindlist(pcas),file.path(out,'PCA.csv'))
fwrite(data.table(component=c('R','edgeR','limma','org.Mm.eg.db'),version=c(as.character(getRversion()),
    vapply(c('edgeR','limma','org.Mm.eg.db'),function(p)as.character(packageVersion(p)),character(1)))),file.path(out,'runtime_R_versions.csv'))
