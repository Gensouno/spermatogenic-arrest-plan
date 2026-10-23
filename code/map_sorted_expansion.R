#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(AnnotationDbi);library(org.Mm.eg.db)})
path<-'data/raw/expansion_assets/GSE290629/GSE290629_RawCountTable_R0032_GEO.tsv.gz'
raw<-read.delim(gzfile(path),check.names=FALSE)
keys<-sub('\\.[0-9]+$','',raw[[1]])
mapping<-mapIds(org.Mm.eg.db,keys=unique(keys),keytype='ENSEMBL',column='SYMBOL',multiVals='asNA')
symbol<-unname(mapping[keys]);keep<-!is.na(symbol)&nzchar(symbol)
x<-as.matrix(raw[-1]);storage.mode(x)<-'double'
stopifnot(!anyDuplicated(colnames(x)),all(is.finite(x)),all(x>=0),all(abs(x-round(x))<1e-6))
counts<-rowsum(x[keep,,drop=FALSE],symbol[keep],reorder=FALSE)
dir.create('data/interim/expansion',recursive=TRUE,showWarnings=FALSE)
con<-gzfile('data/interim/expansion/H3mm15_sorted_symbol_counts.tsv.gz','wt')
write.table(data.frame(gene=rownames(counts),counts,check.names=FALSE),con,sep='\t',quote=FALSE,row.names=FALSE)
close(con)
fwrite(data.frame(input_id=raw[[1]],symbol=symbol,retained=keep),'results/data_expansion/H3mm15_sorted_gene_mapping.csv')
fwrite(data.frame(sample=colnames(x),retained_count_fraction=colSums(counts)/colSums(x)),
       'results/data_expansion/H3mm15_sorted_mapping_mass.csv')
capture.output(sessionInfo(),file='results/data_expansion/sorted_R_sessionInfo.txt')
