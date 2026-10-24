#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)
  library(AnnotationDbi);library(org.Mm.eg.db)})
out<-'results/followup';dir.create(file.path(out,'de'),recursive=TRUE,showWarnings=FALSE)
p<-'data/raw/followup/GSE180204/GSE180204_Lrh1_p1_p3_RnaSeq_Agrimson.rdata.gz'
# Two gzip layers were verified by file signatures. Never execute plot_fpkm or
# reuse the serialized DESeq2 fit/design; extract plain count and metadata slots.
outer<-gzfile(p,'rb');bytes<-readBin(outer,'raw',n=20*1024^2);close(outer)
con<-gzcon(rawConnection(bytes));e<-new.env(parent=emptyenv());objects<-load(con,envir=e);close(con)
stopifnot(setequal(objects,c('dds','plot_fpkm')))
d<-get('dds',e);stopifnot(inherits(d,'DESeqDataSet'))
assays<-attr(attr(attr(d,'assays'),'data'),'listData')
stopifnot(identical(names(assays),'counts'))
x<-assays[['counts']]
m<-as.data.frame(attr(attr(d,'colData'),'listData'))
source_ids<-attr(attr(d,'colData'),'rownames')
stopifnot(is.matrix(x),identical(dim(x),c(55401L,16L)),
  !anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(is.finite(x)),all(x>=0),
  all(abs(x-round(x))<1e-6),all(colSums(x)>0),all(colSums(x)==m$total_counts),
  identical(sub('\\.Aligned.*$','',colnames(x)),source_ids),!anyDuplicated(m$sample_name))
official<-fread(file.path(out,'Lrh1_source_samples.csv'))
stopifnot(setequal(m$sample_name,official$author_sample_name))
official<-official[match(m$sample_name,author_sample_name)]
stopifnot(identical(tolower(official$age),as.character(m$age)),
  identical(ifelse(official$group=='Control','noCre',ifelse(official$age=='P1','NanosCre','DhhCre')),as.character(m$cre)))
meta<-data.frame(official,batch=as.character(m$batch),source_column=colnames(x))
keys<-sub('\\.[0-9]+$','',rownames(x))
mapped<-mapIds(org.Mm.eg.db,keys=unique(keys),keytype='ENSEMBL',column='SYMBOL',multiVals='asNA')
symbols<-unname(mapped[keys]);use<-!is.na(symbols)&nzchar(symbols)
# Preserve every original gene count for filtering/normalization/DE. Annotation
# coverage is diagnostic, not permission to discard unmapped library mass.
counts<-x;rownames(counts)<-keys;stopifnot(!anyDuplicated(keys))
mass<-colSums(x[use,,drop=FALSE])/colSums(x)
fwrite(data.frame(source_column=colnames(x),symbol_mapped_count_fraction=mass),file.path(out,'Lrh1_mapping_mass.csv'))
missing<-data.frame(gene_id=keys[!use],mean_count=rowMeans(x[!use,,drop=FALSE]))
fwrite(missing[order(-missing$mean_count),],file.path(out,'Lrh1_unmapped_genes.csv'))
colnames(counts)<-meta$gsm
fwrite(data.frame(input_id=rownames(x),symbol=symbols,retained=use),file.path(out,'Lrh1_mapping.csv'))
fwrite(data.frame(meta,raw_gene_total=colSums(x),symbol_mapped_count_fraction=mass),file.path(out,'Lrh1_samples.csv'))
fit<-function(counts,meta,adjust=TRUE){
  meta$group<-factor(meta$group,levels=c('Control','KO'));meta$batch<-factor(meta$batch)
  use_batch<-adjust&&nlevels(meta$batch)>1
  design<-if(use_batch) model.matrix(~batch+group,meta) else model.matrix(~group,meta)
  stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
  d<-DGEList(counts=counts);keep<-filterByExpr(d,design=design)
  d<-calcNormFactors(d[keep,,keep.lib.sizes=FALSE])
  v<-voom(d,design,plot=FALSE);f<-eBayes(lmFit(v,design))
  list(table=topTable(f,coef='groupKO',number=Inf,sort.by='none'),expression=v$E,
       design=design,model=if(use_batch)'~ batch + group' else '~ group')
}
summaries<-list();loo<-list();loo_effect<-list();pcas<-list();batch_sensitivity<-list();methods<-list()
for(cohort in unique(meta$cohort)){
  mm<-meta[meta$cohort==cohort,];xx<-counts[,mm$gsm,drop=FALSE]
  stopifnot(sum(mm$group=='KO')==4,sum(mm$group=='Control')==4)
  f<-fit(xx,mm);tt<-f$table
  fwrite(data.frame(gene=rownames(tt),symbol=unname(mapped[rownames(tt)]),tt,row.names=NULL),file.path(out,'de',paste0(cohort,'.csv')))
  # Diagnostic added after observing significance disagreement with the paper;
  # never replaces the prespecified limma primary analysis.
  q<-DGEList(counts=xx[rownames(tt),,drop=FALSE]);q<-calcNormFactors(q)
  q<-estimateDisp(q,f$design,robust=TRUE)
  qfit<-glmQLFit(q,f$design,robust=TRUE)
  qtt<-topTags(glmQLFTest(qfit,coef='groupKO'),n=Inf,sort.by='none')$table
  fwrite(data.frame(gene=rownames(qtt),symbol=unname(mapped[rownames(qtt)]),qtt,row.names=NULL),
         file.path(out,'de',paste0(cohort,'_edgeRQL.csv')))
  methods[[cohort]]<-data.frame(cohort=cohort,primary_method='limma_voom',sensitivity_method='edgeR_robust_QL',
    tested_genes=nrow(tt),limma_DEG=sum(tt$adj.P.Val<.05&abs(tt$logFC)>1),
    edgeRQL_DEG=sum(qtt$FDR<.05&abs(qtt$logFC)>1),
    spearman_logFC=cor(tt[rownames(qtt),'logFC'],qtt$logFC,method='spearman'))
  fwrite(data.frame(gsm=mm$gsm,f$design),file.path(out,paste0(cohort,'_design.csv')))
  summaries[[cohort]]<-data.frame(cohort=cohort,age=mm$age[1],Control=4,KO=4,model=f$model,
    control_genotype=mm$genotype[mm$group=='Control'][1],KO_genotype=mm$genotype[mm$group=='KO'][1],
    tested_genes=nrow(tt),DEG_FDR05_absLFC1=sum(tt$adj.P.Val<.05&abs(tt$logFC)>1),
    min_symbol_mapped_count_fraction=min(mass[meta$cohort==cohort]))
  pc<-prcomp(t(f$expression));vv<-pc$sdev^2/sum(pc$sdev^2)
  pcas[[cohort]]<-data.frame(cohort=cohort,gsm=mm$gsm,group=mm$group,batch=mm$batch,
    pc$x[,1:2],PC1_variance=vv[1],PC2_variance=vv[2])
  for(j in seq_len(ncol(xx))){
    refit<-fit(xx[,-j,drop=FALSE],mm[-j,])$table;common<-intersect(rownames(tt),rownames(refit))
    loo[[paste(cohort,j)]]<-data.frame(cohort=cohort,deleted=mm$gsm[j],genes=length(common),
      spearman=cor(tt[common,'logFC'],refit[common,'logFC'],method='spearman'))
    loo_effect[[paste(cohort,j)]]<-data.frame(cohort=cohort,deleted=mm$gsm[j],gene=rownames(refit),logFC=refit$logFC)
  }
  if(mm$age[1]=='P1'){
    for(mode in c('unadjusted','B2_only')){
      keep<-if(mode=='B2_only')mm$batch=='B2' else rep(TRUE,nrow(mm))
      alt<-fit(xx[,keep,drop=FALSE],mm[keep,],adjust=FALSE)$table
      fwrite(data.frame(gene=rownames(alt),symbol=unname(mapped[rownames(alt)]),alt,row.names=NULL),file.path(out,'de',paste0(cohort,'_',mode,'.csv')))
      common<-intersect(rownames(tt),rownames(alt))
      batch_sensitivity[[mode]]<-data.frame(cohort=cohort,variant=mode,samples=sum(keep),genes=length(common),
        DEG_FDR05_absLFC1=sum(alt$adj.P.Val<.05&abs(alt$logFC)>1),
        spearman_with_batch_adjusted=cor(tt[common,'logFC'],alt[common,'logFC'],method='spearman'))
    }
  }
  print(summaries[[cohort]])
}
fwrite(rbindlist(summaries),file.path(out,'cohort_summary.csv'))
fwrite(rbindlist(loo),file.path(out,'leave_one_out.csv'))
fwrite(rbindlist(loo_effect),file.path(out,'leave_one_out_effects.csv.gz'))
fwrite(rbindlist(pcas),file.path(out,'PCA.csv'))
fwrite(rbindlist(batch_sensitivity),file.path(out,'batch_sensitivity.csv'))
fwrite(rbindlist(methods),file.path(out,'method_sensitivity.csv'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
