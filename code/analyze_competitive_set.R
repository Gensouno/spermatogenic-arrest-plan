#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
out<-commandArgs(trailingOnly=TRUE)[1]
meta<-fread(file.path(out,'prepared','samples.csv'))
annotation<-fread('independent_contexts_snapshot/gene_annotation.csv.gz')
old<-fread('independent_contexts_snapshot/full_statistics.csv.gz')
primary<-fread('replicability_snapshot/all_primary_gene_stage_context.csv.gz')$gene
priorities<-fread('replicability_snapshot/stable_descriptive_priorities.csv')$gene
stopifnot(length(primary)==9675L,!anyDuplicated(primary),length(priorities)==36L,!anyDuplicated(priorities))
results<-list();universes<-list();coverage<-list();audit<-list()
for (cohort_name in unique(meta$cohort)) {
  mm<-meta[cohort==cohort_name]
  d<-read.delim(gzfile(file.path(out,'prepared',paste0(cohort_name,'.tsv.gz'))),check.names=FALSE)
  x<-as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double'
  stopifnot(identical(colnames(x),mm$matrix_sample),ncol(x)==6L,all(x>=0),all(is.finite(x)),all(x==floor(x)))
  group<-factor(mm$group,levels=c('Control','Perturbation'));design<-model.matrix(~group)
  stopifnot(identical(as.integer(table(group)),c(3L,3L)),qr(design)$rank==2)
  y<-DGEList(x);keep<-filterByExpr(y,design=design)
  y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,design,plot=FALSE)
  fit<-eBayes(lmFit(v,design));tt<-topTable(fit,coef=2,number=Inf,sort.by='none')
  previous<-old[cohort==cohort_name]
  stopifnot(identical(rownames(tt),previous$gene))
  errors<-vapply(c('logFC','AveExpr','t','P.Value'),function(k)max(abs(tt[[k]]-previous[[k]])),numeric(1))
  stopifnot(all(errors<1e-8))
  a<-annotation[cohort==cohort_name];a<-a[match(rownames(v),a$gene)]
  stopifnot(identical(a$gene,rownames(v)))
  eligible<-a$globally_unique_symbol & !is.na(a$symbol) & a$symbol%in%primary
  vv<-v[eligible,];aa<-a[eligible];stopifnot(!anyDuplicated(aa$symbol))
  iset<-which(aa$symbol%in%priorities);stopifnot(length(iset)>1L,length(iset)<nrow(vv))
  universes[[cohort_name]]<-data.table(cohort=cohort_name,aa,priority=aa$symbol%in%priorities,
    logFC=tt$logFC[eligible],AveExpr=tt$AveExpr[eligible])
  coverage[[cohort_name]]<-data.table(cohort=cohort_name,symbol=priorities,
    available=priorities%in%aa$symbol,
    status=ifelse(priorities%in%aa$symbol,'tested_unique_primary_symbol','not_tested_or_not_globally_unique'))
  audit[[cohort_name]]<-data.table(cohort=cohort_name,source_tested_genes=nrow(v),universe_genes=nrow(vv),
    measured_priorities=length(iset),background_genes=nrow(vv)-length(iset),libraries=6L,controls=3L,perturbations=3L,
    residual_df=4L,logFC_max_error=errors['logFC'],AveExpr_max_error=errors['AveExpr'],
    t_max_error=errors['t'],PValue_max_error=errors['P.Value'])
  for (policy in c('estimated','fixed_0p01')) {
    rho<-if(policy=='estimated') NA_real_ else .01
    z<-camera(vv,list(priorities=iset),design=design,contrast=2,inter.gene.cor=rho,
      allow.neg.cor=FALSE,use.ranks=FALSE,trend.var=FALSE,directional=TRUE,sort=FALSE)
    reverse_design<-design;reverse_design[,2]<- -reverse_design[,2]
    rev<-camera(vv,list(priorities=iset),design=reverse_design,contrast=2,inter.gene.cor=rho,
      allow.neg.cor=FALSE,use.ranks=FALSE,trend.var=FALSE,directional=TRUE,sort=FALSE)
    stopifnot(abs(z$PValue-rev$PValue)<1e-10,z$Direction!=rev$Direction)
    correlation<-if(policy=='estimated')z$Correlation else rho
    used_vif<-max(1,1+(length(iset)-1)*correlation)
    results[[paste(cohort_name,policy)]]<-data.table(cohort=cohort_name,policy=policy,
      fixed_priorities=36L,NGenes=z$NGenes,universe_genes=nrow(vv),background_genes=nrow(vv)-length(iset),
      correlation_raw_or_assumed=correlation,variance_inflation_used=used_vif,
      camera_df=if(policy=='estimated')min(4L,nrow(vv)-2L) else nrow(vv)-2L,
      Direction=z$Direction,PValue=z$PValue,reverse_PValue=rev$PValue,reverse_Direction=rev$Direction,
      priority_mean_logFC=mean(tt$logFC[eligible][iset]),background_mean_logFC=mean(tt$logFC[eligible][-iset]))
  }
}
z<-rbindlist(results);stopifnot(nrow(z)==8L)
z[,BY:=p.adjust(PValue,'BY')];z[,BH:=p.adjust(PValue,'BH')]
fwrite(z,file.path(out,'camera_results.csv'))
fwrite(rbindlist(universes),file.path(out,'gene_universe.csv.gz'),compress='gzip')
fwrite(rbindlist(coverage),file.path(out,'priority_coverage.csv'))
fwrite(rbindlist(audit),file.path(out,'fit_reconciliation.csv'))
fwrite(data.table(component=c('R','edgeR','limma','data.table'),version=c(as.character(getRversion()),
  vapply(c('edgeR','limma','data.table'),function(p)as.character(packageVersion(p)),character(1)))),file.path(out,'runtime_R_versions.csv'))
writeLines(deparse(getS3method('camera','default',envir=asNamespace('limma'))),file.path(out,'installed_camera_source.txt'))
