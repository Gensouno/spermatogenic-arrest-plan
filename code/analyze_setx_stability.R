#!/usr/bin/env Rscript
# Same-source exploratory count omissions and fixed-universe diagnostics.
source('code/analyze_setx_length.R')

matched_count_pair<-function(counts,lengths,meta) {
  stopifnot(identical(dimnames(counts),dimnames(lengths)),identical(colnames(counts),meta$gsm),
    all(is.finite(counts)),all(counts>=0),all(is.finite(lengths)),all(lengths>0))
  design<-make_design(meta);contrast<-make_contrasts()
  normalized<-normLibSizes(DGEList(counts),method='TMM')
  adjusted<-DGEListFromTximport(list(counts=counts,length=lengths,countsFromAbundance='no'))
  adjusted$samples<-normalized$samples
  zero<-adjusted;zero$offset.prior[]<-0
  stopifnot(identical(adjusted$counts,zero$counts),identical(adjusted$samples,zero$samples))
  libpart<-getOffset(adjusted)-adjusted$offset.prior
  offset_error<-max(abs(libpart-getOffset(zero)))
  stopifnot(offset_error<1e-12)
  result<-list()
  for(name in c('matched_count_length_offset','matched_count_zero_offset_counterfactual')) {
    y<-if(name=='matched_count_length_offset')adjusted else zero
    fit<-glmQLFit(y,design,robust=TRUE,legacy=FALSE,prior.count=.125)
    test<-glmQLFTest(fit,contrast=contrast$interaction)
    stopifnot(all(test$df.test==4))
    delta<-fit$coefficients %*% contrast$interaction/log(2)
    d<-data.table(gene=rownames(test$table),F=test$table$F,P=test$table$PValue,
      BY=p.adjust(test$table$PValue,'BY'),BH=p.adjust(test$table$PValue,'BH'),delta_P18=delta[,'P18_vs_P10'])
    result[[name]]<-d
  }
  list(tables=result,normalization_error=offset_error,
    factors=data.table(gsm=meta$gsm,library_size=normalized$samples$lib.size,TMM_factor=normalized$samples$norm.factors))
}

matched_TPM<-function(tpm,meta) {
  stopifnot(identical(colnames(tpm),meta$gsm),all(is.finite(tpm)),all(tpm>=0))
  design<-make_design(meta);contrasts<-make_contrasts()
  fit<-eBayes(contrasts.fit(lmFit(log2(tpm+.5),design),contrasts$interaction),trend=TRUE,robust=TRUE)
  data.table(gene=rownames(tpm),F=as.numeric(fit$F),P=as.numeric(fit$F.p.value),
    BY=p.adjust(fit$F.p.value,'BY'),BH=p.adjust(fit$F.p.value,'BH'),delta_P18=fit$coefficients[,'P18_vs_P10'])
}

main<-function() {
  out<-'results/setx_stability';meta<-fread('setx_interaction_snapshot/sample_design.csv')
  read<-function(name){d<-read.delim(gzfile(file.path(out,paste0(name,'.tsv.gz'))),check.names=FALSE);x<-as.matrix(d[-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double';x}
  counts<-read('expected_count');lengths<-read('effective_length');tpm<-read('TPM')
  stopifnot(identical(dimnames(counts),dimnames(tpm)),identical(dimnames(counts),dimnames(lengths)))
  map<-fread('new_models_snapshot/gene_annotation.csv')[cohort_source=='Setx']
  stopifnot(identical(map$gene,rownames(counts)))
  good<-!is.na(map$symbol)&nzchar(map$symbol)&!duplicated(map$symbol)&!duplicated(map$symbol,fromLast=TRUE)
  symbols<-setNames(ifelse(good,map$symbol,NA_character_),map$gene)
  candidate_symbols<-fread('replicability_snapshot/candidate_stage_context.csv')$gene
  variants<-c(list(full=integer(),exclude_pair=which(meta$gsm %in% c('GSM7034763','GSM7034783'))),
    setNames(as.list(seq_len(nrow(meta))),paste0('omit_',meta$gsm)))
  full<-list();candidate_rows<-list();summary<-list();agreement<-list()
  previous<-fread('setx_length_snapshot/all_gene_statistics.csv.gz')
  for(name in names(variants)) {
    keep<-!seq_len(nrow(meta)) %in% variants[[name]]
    f<-length_fit(counts[,keep,drop=FALSE],lengths[,keep,drop=FALSE],meta[keep])
    d<-f$table;d[,`:=`(variant=name,symbol=unname(symbols[gene]))]
    candidate_rows[[name]]<-d[!is.na(symbol)&symbol %in% candidate_symbols]
    full[[name]]<-d[,.(variant,gene,F,P,BY,BH,delta_P18)]
    f$summary[,variant:=name];summary[[name]]<-f$summary
    if(name %in% c('full','exclude_pair')) {
      oldname<-if(name=='full')'length_offset_40' else 'length_offset_exclude_pair_38'
      old<-previous[variant==oldname];stopifnot(identical(d$gene,old$gene))
      error<-max(abs(as.matrix(d[,.(F,P,BY,BH,delta_P18)])-as.matrix(old[,.(F,P,BY,BH,delta_P18)])))
      stopifnot(error<1e-9)
      agreement[[name]]<-data.table(variant=name,genes=nrow(d),maximum_absolute_error=error)
    }
    if(startsWith(name,'omit_'))message('Completed ',name)
  }
  fwrite(rbindlist(full),file.path(out,'all_gene_omissions.csv.gz'),compress='gzip')
  fwrite(rbindlist(candidate_rows),file.path(out,'candidate_omissions.csv.gz'),compress='gzip')
  fwrite(rbindlist(summary),file.path(out,'omission_models.csv'))
  fwrite(rbindlist(agreement),file.path(out,'previous_fit_agreement.csv'))
  old_tpm<-fread('setx_interaction_snapshot/all_gene_statistics.csv.gz')[variant=='main']
  old_count<-previous[variant=='length_offset_40']
  common<-rownames(counts)[rownames(counts) %in% intersect(old_tpm$gene,old_count$gene)]
  stopifnot(length(common)==length(intersect(old_tpm$gene,old_count$gene)))
  pair<-matched_count_pair(counts[common,,drop=FALSE],lengths[common,,drop=FALSE],meta)
  pair$tables[['matched_TPM_limma']]<-matched_TPM(tpm[common,,drop=FALSE],meta)
  diagnostics<-rbindlist(lapply(names(pair$tables),function(name){d<-pair$tables[[name]];d[,`:=`(variant=name,symbol=unname(symbols[gene]))];d}))
  fwrite(diagnostics,file.path(out,'matched_diagnostics.csv.gz'),compress='gzip')
  fwrite(pair$factors,file.path(out,'matched_count_normalization.csv'))
  fwrite(data.table(common_tested_genes=length(common),count_library_normalization_max_error=pair$normalization_error),file.path(out,'matched_design_audit.csv'))
  fwrite(data.table(component=c('R','edgeR','limma'),version=c(as.character(getRversion()),as.character(packageVersion('edgeR')),as.character(packageVersion('limma')))),file.path(out,'runtime_versions.csv'))
}
if(sys.nframe()==0L)main()
