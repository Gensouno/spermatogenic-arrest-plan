#!/usr/bin/env Rscript
# Native edgeR v4 length-offset sensitivity on RSEM estimated counts, never rounded.
suppressPackageStartupMessages({library(edgeR);library(data.table)})
source('code/analyze_setx_interaction.R')
length_fit<-function(counts,lengths,meta) {
  stopifnot(identical(dimnames(counts),dimnames(lengths)),identical(colnames(counts),meta$gsm),
    all(is.finite(counts)),all(counts>=0),all(is.finite(lengths)),all(lengths>=0))
  design<-make_design(meta)
  preliminary<-filterByExpr(DGEList(counts),design=design)
  positive<-rowSums(lengths<=0)==0
  keep<-preliminary&positive
  txi<-list(counts=counts[keep,,drop=FALSE],length=lengths[keep,,drop=FALSE],countsFromAbundance='no')
  y<-DGEListFromTximport(txi)
  # Audit the native offset and estimated-count preservation independently.
  expected<-log(txi$length)-rowMeans(log(txi$length))
  stopifnot(max(abs(expected-y$offset.prior))<1e-12,identical(y$counts,txi$counts))
  offset_error<-max(abs(expected-y$offset.prior))
  retained<-filterByExpr(y,design=design);y<-y[retained,,keep.lib.sizes=FALSE]
  y<-normLibSizes(y,method='TMM')
  fit<-glmQLFit(y,design,robust=TRUE,legacy=FALSE,prior.count=.125)
  contrast<-make_contrasts();test<-glmQLFTest(fit,contrast=contrast$interaction)
  stopifnot(all(test$df.test==4))
  effects<-fit$coefficients %*% contrast$genotype/log(2)
  deltas<-fit$coefficients %*% contrast$interaction/log(2)
  result<-data.table(gene=rownames(test$table),F=test$table$F,P=test$table$PValue,
    BY=p.adjust(test$table$PValue,'BY'),BH=p.adjust(test$table$PValue,'BH'))
  for(a in ages)set(result,j=paste0('effect_',a),value=effects[,a])
  for(a in ages[-1])set(result,j=paste0('delta_',a),value=deltas[,paste0(a,'_vs_P10')])
  list(table=result,summary=data.table(samples=ncol(counts),design_rank=qr(design)$rank,interaction_df=4,
    preliminary_expressed=sum(preliminary),excluded_nonpositive_length=sum(preliminary&!positive),
    second_filter_excluded=sum(!retained),tested_genes=nrow(result),
    conditional_BY05=sum(result$BY<.05),conditional_BH05=sum(result$BH<.05),
    noninteger_source_values=sum(abs(counts-round(counts))>1e-6),
    length_offset_max_error=offset_error,
    counts_rounded=FALSE))
}
main<-function(){
  out<-'results/setx_length';meta<-fread('setx_interaction_snapshot/sample_design.csv')
  read<-function(name){d<-read.delim(gzfile(file.path(out,name)),check.names=FALSE);x<-as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double';x}
  counts<-read('expected_count.tsv.gz');lengths<-read('effective_length.tsv.gz')
  map<-fread('new_models_snapshot/gene_annotation.csv')[cohort_source=='Setx']
  stopifnot(identical(map$gene,rownames(counts)),!anyDuplicated(map$gene))
  good<-!is.na(map$symbol)&nzchar(map$symbol)&!duplicated(map$symbol)&!duplicated(map$symbol,fromLast=TRUE)
  symbols<-setNames(ifelse(good,map$symbol,NA_character_),map$gene)
  tables<-list();summaries<-list()
  for(v in c('length_offset_40','length_offset_exclude_pair_38')){
    keep<-if(v=='length_offset_40')rep(TRUE,nrow(meta)) else !meta$gsm %in% c('GSM7034763','GSM7034783')
    f<-length_fit(counts[,keep,drop=FALSE],lengths[,keep,drop=FALSE],meta[keep])
    f$table[,`:=`(variant=v,symbol=unname(symbols[gene]))];f$summary[,variant:=v]
    tables[[v]]<-f$table;summaries[[v]]<-f$summary
  }
  fwrite(rbindlist(tables),file.path(out,'all_gene_statistics.csv.gz'),compress='gzip')
  fwrite(rbindlist(summaries),file.path(out,'variant_summary.csv'))
  fwrite(data.table(component=c('R','edgeR','limma'),version=c(as.character(getRversion()),as.character(packageVersion('edgeR')),as.character(packageVersion('limma')))),file.path(out,'runtime_versions.csv'))
}
if(sys.nframe()==0L)main()
