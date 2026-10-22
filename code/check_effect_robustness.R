#!/usr/bin/env Rscript
# Descriptive sample/parameter sensitivity; no new independent validation claim.
suppressPackageStartupMessages({library(data.table);library(limma);library(edgeR)})
out<-'results/robustness'
dir.create(out,showWarnings=FALSE,recursive=TRUE)
counts<-c('GSE236388','GSE283904','GSE68283','GSE68284','GSE101511')
cohorts<-c(counts,'GSE294209','GSE38096','GSE42601')
read_matrix<-function(path) {
  d<-read.delim(gzfile(path),check.names=FALSE)
  x<-as.matrix(d[-1]); rownames(x)<-d[[1]]; storage.mode(x)<-'double'; x
}
baselines<-lapply(cohorts,function(c) {
  d<-fread(file.path('results/de',paste0(c,'_de.csv')))
  setNames(d$logFC,d$gene)
}); names(baselines)<-cohorts
common<-Reduce(intersect,lapply(baselines,names))
fit_effect<-function(x,group,kind,pseudocount=.5) {
  group<-factor(group,levels=c('WT','KO')); design<-model.matrix(~group)
  if (kind=='counts') {
    d<-DGEList(counts=x,group=group)
    keep<-filterByExpr(d,design=design)
    d<-calcNormFactors(d[keep,,keep.lib.sizes=FALSE])
    fit<-eBayes(lmFit(voom(d,design,plot=FALSE),design))
  } else {
    if(kind=='TPM') {
      x<-log2(x[rowSums(x>=1)>=min(3,ncol(x)),,drop=FALSE]+pseudocount)
    }
    fit<-eBayes(lmFit(x,design),trend=TRUE,robust=TRUE)
  }
  fit$coefficients[,'groupKO']
}
metrics<-list(); pairwise<-list(); details<-list(); checks<-list()
record<-function(gse,variant,effect,kind,dropped='',min_group=NA_integer_) {
  ref<-baselines[[gse]]; shared<-intersect(names(ref),names(effect))
  strong<-shared[abs(ref[shared])>.5]
  metrics[[length(metrics)+1]]<<-data.frame(gse=gse,variant=variant,kind=kind,dropped=dropped,
    min_group_n=min_group,compared_genes=length(shared),
    spearman_with_full=cor(ref[shared],effect[shared],method='spearman'),
    strong_full_genes=length(strong),
    sign_agreement_strong=if(length(strong)) mean(sign(ref[strong])==sign(effect[strong])) else NA_real_)
  for(other in setdiff(cohorts,gse)) {
    genes<-intersect(common,names(effect))
    reference<-baselines[[other]][genes]
    pairwise[[length(pairwise)+1]]<<-data.frame(gse=gse,other=other,variant=variant,kind=kind,dropped=dropped,
        genes=length(genes),baseline_rho=cor(ref[genes],reference,method='spearman'),
        variant_rho=cor(effect[genes],reference,method='spearman'))
  }
}
for(gse in cohorts) {
  kind<-if(gse %in% counts) 'counts' else if(gse=='GSE294209') 'TPM' else 'array'
  path<-if(kind=='array') file.path('results/arrays',paste0(gse,'_log2_expression.tsv.gz')) else
                            file.path('data/bulk',paste0(gse,'_matrix.tsv.gz'))
  meta_path<-if(kind=='array') file.path('results/arrays',paste0(gse,'_samples.tsv')) else
                            file.path('data/bulk',paste0(gse,'_samples.tsv'))
  x<-read_matrix(path); meta<-fread(meta_path)
  stopifnot(identical(colnames(x),meta$matrix_sample))
  full<-fit_effect(x,meta$group,kind)
  genes<-intersect(names(full),names(baselines[[gse]]))
  err<-max(abs(full[genes]-baselines[[gse]][genes]))
  stopifnot(setequal(names(full),names(baselines[[gse]])),err<1e-8)
  checks[[gse]]<-data.frame(gse=gse,baseline_max_error=err)
  for(i in seq_len(ncol(x))) {
    effect<-fit_effect(x[,-i,drop=FALSE],meta$group[-i],kind)
    variant<-paste0('drop_',meta$matrix_sample[i])
    record(gse,variant,effect,'sample_leave_one_out',meta$matrix_sample[i],min(table(meta$group[-i])))
    details[[length(details)+1]]<-data.table(gse=gse,variant=variant,gene=names(effect),logFC=as.numeric(effect))
  }
  if(kind=='TPM') {
    for(p in c(.1,1,2)) {
      effect<-fit_effect(x,meta$group,kind,pseudocount=p)
      record(gse,paste0('TPM_pseudocount_',p),effect,'TPM_pseudocount')
    }
  }
  if(kind=='array') {
    probe<-readRDS(file.path('results/arrays',paste0(gse,'_normalized_probes.rds')))
    # Same normalized probes, alternative gene summaries. Selection uses mean
    # intensity only, not case/control difference or a significance statistic.
    e<-probe$expression; symbols<-probe$symbols
    by_gene<-split(seq_len(nrow(e)),symbols)
    med<-t(vapply(by_gene,function(ix) apply(e[ix,,drop=FALSE],2,median),numeric(ncol(e))))
    best<-vapply(by_gene,function(ix) ix[which.max(rowMeans(e[ix,,drop=FALSE]))],integer(1))
    representative<-e[best,,drop=FALSE]; rownames(representative)<-names(best)
    for(method in c('median','max_mean_probe')) {
      alternative<-if(method=='median') med else representative
      effect<-fit_effect(alternative,meta$group,'array')
      record(gse,method,effect,'array_probe_summary')
    }
  }
  cat(gse,'sample sensitivity done\n')
}
fwrite(rbindlist(metrics),file.path(out,'sensitivity_metrics.csv'))
fwrite(rbindlist(pairwise),file.path(out,'pairwise_sensitivity.csv'))
fwrite(rbindlist(checks),file.path(out,'baseline_reproduction.csv'))
fwrite(rbindlist(details),file.path(out,'leave_one_out_logFC.csv.gz'),compress='gzip')
writeLines(c('No samples automatically removed.',
  'LOO counts rerun expression filtering, TMM and voom; TPM reruns abundance filtering and limma-trend.',
  'Array LOO conditions on full-cohort normalization; it is not full preprocessing leave-one-out.',
  'Sin3a LOO leaves one group n=1: a stress diagnostic, not reliable inferential replication.',
  'TPM filter unchanged in full-data pseudocount comparisons.',
  'Probe summary sensitivity conditions on fixed normalization.',
  'Pairwise sensitivity varies one study at a time; ranges are not confidence intervals.'),file.path(out,'interpretation.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
