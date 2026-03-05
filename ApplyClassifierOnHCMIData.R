source("CodeBlocks/convertHCMIDataEntrez.R")
source("CodeBlocks/getHCMI_tissue_expression.R")
source("CodeBlocks/RemoveHCMIPurity.R")
source("CodeBlocks/GetBestModelParam.R")
source("CodeBlocks/multiclass_pairs_classifier.R")
source("CodeBlocks/convertHCMIDataEntrez.R")
source("CodeBlocks/getHCMI_tissue_expression.R")
source("CodeBlocks/getTissueExpression_TCGA.R")
source("CodeBlocks/RemovePurity.R")

dir.create("ClassifierResults")
tissues=c("BLCA", "PAAD", "SARC", "OV","COADREAD", "LIHCCHOL","LUAD", "UCEC", "BRCA", 
          "GEA", "LGGGBM","SKCM", "LUSC", "OV", "KIRCKICH" )

hcmi_dat=convertHCMIData()
for (tissue in tissues){ 
  tissue_hcmi=getHCMITissueExp(hcmi_dat, tissue)
  adjusted_tissue_hcmi=HCMI_PurityCorrect(tissue_hcmi)
  tissue_tmp=GetTMPTissueExpression(tissue)
  rel_tissue_tcga=tissue_tmp$gexp[intersect(rownames(tissue_tmp$gexp), rownames(tissue_hcmi)),]
  
  purity_removed_tcga=AdjustData(rel_tissue_tcga)
  purity_removed_tcga_anno=tissue_tmp$anno[match(colnames(purity_removed_tcga), tissue_tmp$anno[,1]), ] 
  best_param=GetBestModelParameters(tissue)
  
  tm=as.integer(str_split(best_param$TrainA_TestA, "_")[[1]]); ta_ta_x=tm[1];ta_ta_y=tm[2]
  tm1=as.integer(str_split(best_param$TrainNA_TestNA, "_")[[1]]); tna_tna_x=tm1[1];tna_tna_y=tm1[2]
  
  res_ta_ta=MultiClassClassifier(purity_removed_tcga, unlist(purity_removed_tcga_anno[,2]), adjusted_tissue_hcmi, ta_ta_x, ta_ta_y) 
  res_tna_tna=MultiClassClassifier(rel_tissue_tcga,unlist(tissue_tmp$anno[,2]), tissue_hcmi, tna_tna_x, tna_tna_y) 
  saveRDS(file=paste("ClassifierResults/", tissue, "_TA_TA_res.RDS", sep=""),res_ta_ta)
  saveRDS(file=paste("ClassifierResults/", tissue, "_TNA_TNA_res.RDS", sep=""),res_tna_tna)
}



  