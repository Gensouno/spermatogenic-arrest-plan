#!/usr/bin/env Rscript
# Independent R p.adjust check on the entire family, not selected directions.
suppressPackageStartupMessages(library(data.table))
out<-'results/replicability'
d<-fread(file.path(out,'directional_replicability_tests.csv.gz'))
audit<-list()
for(u in unique(d$universe)) {
  part<-d[universe==u]
  stopifnot(nrow(part)==4L*uniqueN(part$gene),!anyDuplicated(part[,.(gene,direction,r)]))
  by_error<-max(abs(p.adjust(part$partial_conjunction_P,'BY')-part$joint_BY_q))
  bh_error<-max(abs(p.adjust(part$partial_conjunction_P,'BH')-part$joint_BH_q_sensitivity))
  stopifnot(is.finite(by_error),is.finite(bh_error),max(by_error,bh_error)<1e-12)
  audit[[u]]<-data.table(universe=u,hypotheses=nrow(part),BY_max_error=by_error,BH_max_error=bh_error)
}
fwrite(rbindlist(audit),file.path(out,'R_adjustment_audit.csv'))
