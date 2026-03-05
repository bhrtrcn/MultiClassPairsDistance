getHCMITissueExp=function(data, tissue){
  library(readr)
  library(readxl)
  
i f(tissue=="GEA"){
    tissue="STADESCA"
  }

  hcmi_cli=read_excel("~/Desktop/HCMI_very_final/HCMI_637_TCGA_CCLE_ClinicalTeam_20240618.xlsx")
  cid=unlist(subset(hcmi_cli, Diagnosis_TMP_Code==tissue, select="Case_ID"))
  tissue_fpkm_uq=data[,sapply(colnames(data), substring, 1,17)%in%cid]
  rm=which(apply(tissue_fpkm_uq, 1,sum)==0)
  if(length(rm)>0){ 
      hcmi_t_gexp=tissue_fpkm_uq[-rm, ]
  }
  return(hcmi_t_gexp)
}
