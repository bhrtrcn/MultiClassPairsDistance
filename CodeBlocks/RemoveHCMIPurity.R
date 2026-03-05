HCMI_PurityCorrect=function(dat){ 
  library(limma)
  purity=read.csv("4.26.24_final_consensus_pp.txt", sep="\t")
  purity_hcmi=purity[, c("sample_id", "consensus_purity")]
  purity_hcmi=purity_hcmi[!duplicated(purity_hcmi$sample_id),]
  incs_=intersect(purity_hcmi$sample_id, colnames(dat))
  rownames(purity_hcmi)=purity_hcmi$sample_id
  dat=log2(dat[,incs_] +1)
  purity_hcmi=purity_hcmi[incs_,]

  tumors_purity=as.data.frame(purity_hcmi)
  colnames(tumors_purity)=c("sample_id", "purity")
  tumors_purity$purity=as.numeric(tumors_purity$purity)
  tumors_infiltrate = 1 - tumors_purity$purity
  design = model.matrix(~ tumors_infiltrate)
  fit <- lmFit(dat,design)
  beta <- fit$coefficients[,2,drop=FALSE]
  dat_adjusted = as.matrix(dat) - beta %*% t(tumors_infiltrate)
  return(dat_adjusted)
}
