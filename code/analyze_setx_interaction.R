#!/usr/bin/env Rscript
# Exploratory conditional source-library model; no animal-level certification.
suppressPackageStartupMessages({library(data.table);library(limma)})
ages <- c('P10','P12','P14','P16','P18')
cells <- as.vector(t(outer(ages,c('Control','KO'),paste,sep='_')))
make_design <- function(meta) {
  cell <- factor(paste(meta$age,meta$group,sep='_'),levels=cells)
  if(anyNA(cell) || any(table(cell)<2)) stop('Invalid or insufficient design cells')
  design <- model.matrix(~0+cell); colnames(design) <- cells
  stopifnot(qr(design)$rank==10)
  design
}
make_contrasts <- function() {
  genotype <- matrix(0,10,5,dimnames=list(cells,ages))
  for(a in ages) {genotype[paste0(a,'_KO'),a]<-1;genotype[paste0(a,'_Control'),a]<- -1}
  interaction <- sapply(ages[-1],function(a)genotype[,a]-genotype[,'P10'])
  colnames(interaction)<-paste0(ages[-1],'_vs_P10')
  stopifnot(qr(interaction)$rank==4)
  list(genotype=genotype,interaction=interaction)
}
expressed <- function(x,meta) {
  cell <- paste(meta$age,meta$group,sep='_')
  Reduce(`|`,lapply(cells,function(k)rowSums(x[,cell==k,drop=FALSE]>=1)>=min(4,sum(cell==k))))
}
fit_interaction <- function(x,meta,pseudo=.5,moderate=TRUE) {
  stopifnot(identical(colnames(x),meta$gsm),all(is.finite(x)),all(x>=0),pseudo>0)
  design<-make_design(meta); keep<-expressed(x,meta); expression<-log2(x[keep,,drop=FALSE]+pseudo)
  fit<-lmFit(expression,design); contrasts<-make_contrasts()
  g<-contrasts.fit(fit,contrasts$genotype)$coefficients
  f<-contrasts.fit(fit,contrasts$interaction)
  # Independently check coefficients against direct cell mean differences.
  means<-sapply(cells,function(k)rowMeans(expression[,paste(meta$age,meta$group,sep='_')==k,drop=FALSE]))
  direct<-means %*% contrasts$interaction
  stopifnot(max(abs(direct-f$coefficients))<1e-10)
  result<-data.table(gene=rownames(expression))
  for(a in ages)set(result,j=paste0('effect_',a),value=g[,a])
  for(a in ages[-1])set(result,j=paste0('delta_',a),value=f$coefficients[,paste0(a,'_vs_P10')])
  f_error<-NA_real_
  if(moderate) {
    f<-eBayes(f,trend=TRUE,robust=TRUE)
    # Independent nested-model residual calculation of the same four-df F.
    reduced_design<-model.matrix(~factor(meta$age)+factor(meta$group))
    stopifnot(qr(reduced_design)$rank==6)
    reduced<-lmFit(expression,reduced_design)
    rss_full<-rowSums((expression-t(design %*% t(fit$coefficients)))^2)
    rss_reduced<-rowSums((expression-t(reduced_design %*% t(reduced$coefficients)))^2)
    independent_F<-pmax(0,(rss_reduced-rss_full)/4)/f$s2.post
    f_error<-max(abs(independent_F-f$F)/pmax(1,abs(f$F)))
    stopifnot(f_error<1e-8)
    result[,`:=`(F=as.numeric(f$F),P=as.numeric(f$F.p.value),
        BY=p.adjust(f$F.p.value,'BY'),BH=p.adjust(f$F.p.value,'BH'))]
  }
  list(table=result,samples=ncol(x),rank=qr(design)$rank,residual_df=ncol(x)-qr(design)$rank,
       direct_coefficient_max_error=max(abs(direct-f$coefficients)),nested_F_max_relative_error=f_error)
}
main <- function() {
  out<-'results/setx_interaction';meta<-fread(file.path(out,'sample_design.csv'))
  d<-read.delim(gzfile(file.path(out,'TPM.tsv.gz')),check.names=FALSE)
  x<-as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double'
  map<-fread('new_models_snapshot/gene_annotation.csv')[cohort_source=='Setx']
  stopifnot(!anyDuplicated(map$gene),identical(map$gene,rownames(x)))
  valid<-!is.na(map$symbol)&nzchar(map$symbol)
  unique_symbol<-valid&!duplicated(map$symbol)&!duplicated(map$symbol,fromLast=TRUE)
  symbol_map<-setNames(ifelse(unique_symbol,map$symbol,NA_character_),map$gene)
  candidates<-fread('replicability_snapshot/candidate_stage_context.csv')$gene
  variants<-list(main=list(pseudo=.5,drop=character()),pseudo_0p1=list(pseudo=.1,drop=character()),
    pseudo_1=list(pseudo=1,drop=character()),exclude_ambiguous_pair=list(pseudo=.5,drop=c('GSM7034763','GSM7034783')))
  summaries<-list();full<-list();candidate_effects<-list()
  for(name in names(variants)) {
    v<-variants[[name]];keep<-!meta$gsm %in% v$drop
    f<-fit_interaction(x[,keep,drop=FALSE],meta[keep],v$pseudo)
    f$table[,symbol:=unname(symbol_map[gene])]
    f$table[,variant:=name];full[[name]]<-f$table
    candidate_effects[[name]]<-f$table[!is.na(symbol)&symbol %in% candidates]
    summaries[[name]]<-data.table(variant=name,pseudocount=v$pseudo,samples=f$samples,design_rank=f$rank,
      residual_df=f$residual_df,interaction_df=4,tested_genes=nrow(f$table),
      conditional_BY05=sum(f$table$BY<.05),conditional_BH05=sum(f$table$BH<.05),
      direct_coefficient_max_error=f$direct_coefficient_max_error,nested_F_max_relative_error=f$nested_F_max_relative_error)
  }
  # Omission fits are effect checks only, with their own expression filters.
  for(i in seq_len(ncol(x))) {
    f<-fit_interaction(x[,-i,drop=FALSE],meta[-i],.5,FALSE)$table
    f[,symbol:=unname(symbol_map[gene])];f[,variant:=paste0('omit_',meta$gsm[i])]
    candidate_effects[[paste0('omit_',i)]]<-f[!is.na(symbol)&symbol %in% candidates]
  }
  fwrite(rbindlist(full),file.path(out,'all_gene_statistics.csv.gz'),compress='gzip')
  fwrite(rbindlist(candidate_effects,fill=TRUE),file.path(out,'candidate_effects.csv.gz'),compress='gzip')
  fwrite(rbindlist(summaries),file.path(out,'variant_summary.csv'))
  fwrite(data.table(component=c('R','limma','data.table','statmod'),version=c(as.character(getRversion()),
    vapply(c('limma','data.table','statmod'),function(p)as.character(packageVersion(p)),character(1)))),
    file.path(out,'runtime_R_versions.csv'))
}
if(sys.nframe()==0L)main()
