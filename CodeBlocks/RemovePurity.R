
AdjustData=function(dat){ 
  library(edgeR)
  dat=rel_tissue_tcga
  for_purity=read.csv2("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
  rels=for_purity[,c("array","purity")]
  incs=which(sapply(rels$array, substring, 13, 16)=="-01")
  rels=rels[incs,]
  rels$sample_id=sapply(rels$array, substring, 1, 12)
  inc_samples=rels[match(colnames(dat), rels$sample_id),c("sample_id", "purity")]
  inc_samples=inc_samples[inc_samples$purity!="",]
  inc_samples=inc_samples[complete.cases(inc_samples),]
  rownames(inc_samples)=inc_samples$sample_id
  incs_=intersect(inc_samples$sample_id, colnames(dat))
  dat=log2(dat[,incs_] +1)
  inc_samples=inc_samples[incs_,]
  tumors_purity=as.data.frame(inc_samples)
  tumors_purity$purity=as.numeric(tumors_purity$purity)
  tumors_infiltrate = 1 - tumors_purity$purity
  design = model.matrix(~ tumors_infiltrate)
  fit <- lmFit(dat,design)
  beta <- fit$coefficients[,2,drop=FALSE]
  dat_adjusted = as.matrix(dat) - beta %*% t(tumors_infiltrate)
  return(dat_adjusted)
}
