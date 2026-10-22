#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(limma); library(affy); library(mouse4302.db); library(mouse4302cdf)
  library(AnnotationDbi); library(data.table)
})
dir.create('results/arrays',showWarnings=FALSE,recursive=TRUE)
audits<-list()
for (gse in c('GSE38096','GSE42601')) {
  base<-file.path('data/interim/arrays',gse)
  meta<-read.delim(file.path(base,'samples.tsv'),check.names=FALSE)
  files<-file.path(base,meta$raw_file)
  if (gse=='GSE38096') {
    # Raw scan metadata and GEO both identify a SINGLE Cy3 channel.
    obj<-read.maimages(files,source='agilent',green.only=TRUE,
                       columns=list(G='gMedianSignal',Gb='gBGMedianSignal'))
    obj<-backgroundCorrect(obj,method='normexp',offset=50)
    norm<-normalizeBetweenArrays(obj,method='quantile')
    e<-norm$E
    probes<-obj$genes$ProbeName
    annotation<-read.delim(file.path(base,'probe_annotation.tsv'),check.names=FALSE)
    annotation<-annotation[!is.na(annotation$NAME)&nzchar(annotation$NAME)&
                            !is.na(annotation$GENE_SYMBOL)&nzchar(annotation$GENE_SYMBOL),]
    mapping<-tapply(annotation$GENE_SYMBOL,annotation$NAME,function(z) {
      u<-unique(z); if(length(u)==1) u else NA_character_
    })
    symbols<-unname(mapping[probes])
    keep<-obj$genes$ControlType==0 & !is.na(symbols) & nzchar(symbols)
    method<-'Agilent single Cy3; normexp offset50; quantile; mean log2 probes per symbol'
  } else {
    raw<-ReadAffy(filenames=files)
    if (annotation(raw)!='mouse4302') stop('Unexpected CEL chip type')
    obj<-rma(raw,verbose=FALSE)
    e<-exprs(obj); probes<-rownames(e)
    mapping<-mapIds(mouse4302.db,keys=probes,keytype='PROBEID',column='SYMBOL',multiVals='asNA')
    symbols<-unname(mapping[probes])
    keep<-!is.na(symbols)&nzchar(symbols)
    method<-'Affymetrix mouse4302; RMA; mean log2 probesets per symbol'
  }
  stopifnot(all(is.finite(e)),length(symbols)==nrow(e))
  fwrite(data.frame(probe=probes,symbol=symbols,retained=keep),
         file.path('results/import_qc',paste0(gse,'_probe_mapping.csv')))
  saveRDS(list(expression=e[keep,,drop=FALSE],symbols=symbols[keep],probes=probes[keep],meta=meta),
          file.path('results/arrays',paste0(gse,'_normalized_probes.rds')))
  expr<-avereps(e[keep,,drop=FALSE],ID=symbols[keep])
  colnames(expr)<-meta$matrix_sample
  # Arrays remain in log2 intensity units; do not feed these into CPM NNLS.
  con<-gzfile(file.path('results/arrays',paste0(gse,'_log2_expression.tsv.gz')),'wt')
  write.table(data.frame(gene=rownames(expr),expr,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE)
  close(con)
  group<-factor(meta$group,levels=c('WT','KO'))
  design<-model.matrix(~group)
  fit<-eBayes(lmFit(expr,design),trend=TRUE,robust=TRUE)
  tt<-topTable(fit,coef='groupKO',number=Inf,sort.by='none')
  fwrite(data.frame(gene=rownames(tt),tt,row.names=NULL),file.path('results/de',paste0(gse,'_de.csv')))
  pca<-prcomp(t(expr),center=TRUE,scale.=FALSE)
  frac<-pca$sdev^2/sum(pca$sdev^2)
  fwrite(data.frame(sample=colnames(expr),group=group,pca$x[,1:2],PC1_variance=frac[1],PC2_variance=frac[2]),
         file.path('results/import_qc',paste0(gse,'_pca.csv')))
  fwrite(meta,file.path('results/arrays',paste0(gse,'_samples.tsv')),sep='\t')
  audits[[gse]]<-data.frame(gse=gse,WT=sum(group=='WT'),KO=sum(group=='KO'),
              input_probes=nrow(e),retained_probes=sum(keep),tested_genes=nrow(expr),
              DEG_FDR05_absLFC1=sum(tt$adj.P.Val<.05&abs(tt$logFC)>1),method=method)
  print(audits[[gse]])
}
fwrite(rbindlist(audits),'results/import_qc/microarray_cohorts.csv')
capture.output(sessionInfo(),file='results/import_qc/microarray_R_sessionInfo.txt')
