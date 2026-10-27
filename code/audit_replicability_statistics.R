#!/usr/bin/env Rscript
# Recompute moderated-t tails from prepared expression, without changing old results.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
out <- 'results/replicability'; dir.create(out,recursive=TRUE,showWarnings=FALSE)
primary <- c('GSE236388','GSE283904','GSE68283','GSE68284','GSE101511')
cohorts <- c(primary,'GSE294209','GSE38096','GSE42601')
records <- list(); audits <- list()
for (gse in cohorts) {
  kind <- if(gse %in% primary) 'counts' else if(gse=='GSE294209') 'TPM' else 'array'
  base <- if(kind=='array') 'results/arrays' else 'data/bulk'
  expression_path <- file.path(base,paste0(gse,if(kind=='array') '_log2_expression.tsv.gz' else '_matrix.tsv.gz'))
  d <- read.delim(gzfile(expression_path),check.names=FALSE)
  x <- as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double'
  meta <- fread(file.path(base,paste0(gse,'_samples.tsv')))
  stopifnot(identical(colnames(x),meta$matrix_sample),!anyDuplicated(rownames(x)),all(is.finite(x)),all(table(meta$group)>=2))
  group<-factor(meta$group,levels=c('WT','KO'));design<-model.matrix(~group)
  if(kind=='counts') {
    stopifnot(all(x>=0),all(abs(x-round(x))<1e-6))
    y<-DGEList(counts=x,group=group);keep<-filterByExpr(y,design=design)
    y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
    fit<-eBayes(lmFit(voom(y,design,plot=FALSE),design))
  } else {
    if(kind=='TPM') x<-log2(x[rowSums(x>=1)>=3,,drop=FALSE]+.5)
    fit<-eBayes(lmFit(x,design),trend=TRUE,robust=TRUE)
  }
  tt<-topTable(fit,coef='groupKO',number=Inf,sort.by='none')
  previous<-fread(file.path('results/de',paste0(gse,'_de.csv')))
  stopifnot(setequal(rownames(tt),previous$gene))
  previous<-previous[match(rownames(tt),gene)]
  errors<-vapply(c('logFC','t','P.Value','adj.P.Val'),function(k) max(abs(tt[[k]]-previous[[k]])),numeric(1))
  stopifnot(all(is.finite(errors)),max(errors)<1e-8)
  t<-fit$t[,'groupKO'];df<-rep_len(fit$df.total,length(t))
  p_up<-pt(t,df=df,lower.tail=FALSE);p_down<-pt(t,df=df,lower.tail=TRUE)
  stopifnot(all(is.finite(df)),all(df>0),all(is.finite(p_up)),all(is.finite(p_down)),
    max(abs(p_up+p_down-1))<1e-12,max(abs(2*pmin(p_up,p_down)-tt$P.Value))<1e-12,
    all(sign(tt$logFC)==sign(t)))
  records[[gse]]<-data.table(gse=gse,gene=rownames(tt),unit=kind,logFC=tt$logFC,t=as.numeric(t),
      df_total=df,P_two_sided=tt$P.Value,original_adj_P=tt$adj.P.Val,P_up=p_up,P_down=p_down)
  audits[[gse]]<-data.table(gse=gse,unit=kind,samples=ncol(x),WT=sum(group=='WT'),KO=sum(group=='KO'),
      tested_genes=nrow(tt),max_logFC_error=errors[['logFC']],max_t_error=errors[['t']],
      max_P_error=errors[['P.Value']],max_adjusted_P_error=errors[['adj.P.Val']],
      min_df=min(df),max_df=max(df),scope='Prepared expression re-fit, not raw-array preprocessing or annotation regeneration')
  cat('Reconstructed source statistics:',gse,nrow(tt),'genes\n')
}
fwrite(rbindlist(records),file.path(out,'recomputed_study_statistics.csv.gz'),compress='gzip')
fwrite(rbindlist(audits),file.path(out,'source_statistics_audit.csv'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
fwrite(data.table(component=c('R','edgeR','limma','data.table'),
  version=c(as.character(getRversion()),vapply(c('edgeR','limma','data.table'),function(p) as.character(packageVersion(p)),character(1)))),
  file.path(out,'runtime_R_versions.csv'))
