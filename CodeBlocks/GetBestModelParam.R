
GetBestModelParameters=function(tissue){
  best_mods=read.csv("model_accuracies.csv")
  tissue_best=subset(best_mods, Tissue==tissue)
  traina_testa=tissue_best[which.max(tissue_best[,"TrainA_TestA"]), "Model"]
  trainna_testna=tissue_best[which.max(tissue_best[,"TrainNA_TestNA"]), "Model"]
  return(list(TrainA_TestA=traina_testa,TrainNA_TestNA= trainna_testna))
}