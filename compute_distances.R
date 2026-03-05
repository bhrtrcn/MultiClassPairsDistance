library(readxl)
library(stringr)
library(reshape2)

`%notin%`=Negate(`%in%`)
dat2=read_excel("Analysis Tracker - Consensus Samples for HCMI.xlsx", "RNA_integrated_distances" )

#exclusion list
dat2=dat2[,c(1,2,3,9,52)]
dat2$Excluded="No"
exc_models=read_excel("excluded_models.xlsx")

# TMP tissue information 
cli=read_excel("HCMI_637_TCGA_CCLE_ClinicalTeam_20240618.xlsx")
tmp_codes=cli[match(dat2$case_id3, cli$Case_ID),"Diagnosis_TMP_Code"]
tmp_codes[tmp_codes=="STADESCA"]="GEA"
dat2$Diagnosis_TMP_Code=unlist(tmp_codes)

dat2[dat2$aliquot_id3%in%exc_models$model_id, "Excluded"]="Yes"

tissues=unique(sapply(list.files("ClassifierResults"), function(x){str_split(x, "_")[[1]][1]}))
ta_distances=list()

for (tissue in tissues){ 
  ta_ta1=readRDS(paste("ClassifierResults/", tissue, "_TA_TA_res.RDS", sep=""))
  ta_ta=ta_ta1$Probabilities
  ta_distance=matrix(nrow=nrow(ta_ta), ncol=nrow(ta_ta))
  for (i in 1:nrow(ta_ta)){
    for (j in 1:nrow(ta_ta)){
      A=ta_ta[i,]
      B=ta_ta[j,]
      ta_distance[i,j] <- sqrt(sum((A - B)^2))
    }
  }
  rownames(ta_distance)=rownames(ta_ta)
  colnames(ta_distance)=rownames(ta_ta)
}

ta_all=do.call(rbind,ta_distances)
colnames(ta_all)=c("Sample1", "Sample2", "ISB_Distance")

merged1=merge(dat2, ta_all,  by.x=c("aliquot_id3","matched_tumor_aliquot"), by.y=c("Sample1", "Sample2"), all.x=TRUE)
merged1[merged1$Excluded=="Yes", "ISB_Distance"]=NA
merged1$ISB_Outlier=NA

all_outliers=c()
for (tissue in tissues){
  t_d=merged1[merged1$Diagnosis_TMP_Code==tissue,c("aliquot_id3", "ISB_Distance")]
  Q1 <- quantile(na.omit(t_d$ISB_Distance), 0.25)
  Q3 <- quantile(na.omit(t_d$ISB_Distance), 0.75)
  IQR <- Q3 - Q1
  upper_bound <- Q3 + 1.5 * IQR
  outliers <-t_d$aliquot_id3[which(t_d$ISB_Distance> upper_bound)]
  all_outliers=c(all_outliers, outliers)
}
merged1[!is.na(merged1$ISB_Distance),"ISB_Outlier"]="No"
merged1[merged1$aliquot_id3%in%all_outliers, "ISB_Outlier"]="Yes"

ff=merged1[match(dat2$aliquot_id3, merged1$aliquot_id3),]
ff$ISB_Outlier[ff$ISB_Outlier=="Yes"]=TRUE
ff$ISB_Outlier[ff$ISB_Outlier=="No"]=FALSE

write.csv(ff,"multiclasspair_distance.csv")
