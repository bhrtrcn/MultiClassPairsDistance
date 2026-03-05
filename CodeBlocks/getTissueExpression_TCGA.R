GetTMPTissueExpression=function(tissue){
  library(stringr)
  library(readr)
  tmp=as.data.frame(read_tsv(paste("TMP_20230209/", tissue, "_v12_20210228.tsv", sep="")))
  dat=tmp[,3:ncol(tmp)]
  gexp=dat[,colnames(dat)[sapply(colnames(dat), startsWith, "N:GEXP")]]
  colnames(gexp)=paste("T", sapply(colnames(gexp), function(x){str_split(x, "\\:")[[1]][5]}), sep="")
  rownames(gexp)=tmp[,tissue]
  gexp[gexp<0]=0
  rmr=which(apply(gexp,1,sum)==0)
  if(length(rmr)>0){
    gexp=gexp[-rmr,]
  }
  return(list(gexp=t(gexp), anno=tmp[,1:2]))
}
  
  