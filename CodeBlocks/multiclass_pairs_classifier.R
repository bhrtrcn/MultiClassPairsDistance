library(multiclassPairs)
MultiClassClassifier=function(train_data, train_labels, test_data, x, y){
  object <- ReadData(Data =train_data, 
                     Labels = as.character(train_labels), 
                     Platform = NULL, 
                     verbose = FALSE) 
  
  genes_RF <- sort_genes_RF(data_object = object,
                            platform_wise = FALSE,
                            num.trees = 5000, 
                            verbose = FALSE)
  
  rules_RF <- sort_rules_RF(data_object = object, 
                            sorted_genes_RF = genes_RF,
                            genes_altogether = x,
                            genes_one_vs_rest = y, 
                            num.trees = 5000,
                            verbose = TRUE)
  
  
  RF_classifier <- train_RF(data_object = object,
                            sorted_rules_RF = rules_RF,
                            rules_altogether = x,
                            rules_one_vs_rest = y,
                            run_boruta = TRUE, 
                            plot_boruta = FALSE,
                            probability = TRUE,
                            num.trees = 5000,
                            boruta_args = list(),
                            verbose = TRUE)
  
  results <- predict_RF(classifier = RF_classifier, 
                        Data =test_data,
                        impute = TRUE)
  
  test_pred <- results$predictions
  if (is.factor(test_pred)) {
    x <- as.character(test_pred)
  }
  # if the classifier trained using probability   = TRUE
  if (is.matrix(test_pred)) {
    x <- colnames(test_pred)[max.col(test_pred)]
  }
  
  return(list(Assignments=cbind.data.frame(Samples=colnames(test_data), Assignment=x), Probabilities=test_pred, Rules= RF_classifier$RF_scheme$rules))
}
