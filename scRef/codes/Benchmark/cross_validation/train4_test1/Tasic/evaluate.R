evaluate <- function(TrueLabelsPath, PredLabelsPath){
  "
  Script to evaluate the performance of the classifier.
  It returns multiple evaluation measures: the confusion matrix, median F1-score, F1-score for each class, accuracy, percentage of unlabeled, population size. 

  The percentage of unlabeled cells is find by checking for cells that are labeled 'Unassigned', 'unassigned', 'Unknown', 'Nodexx', 'rand', or 'ambiguous'.
  
  Parameters
  ----------
  TrueLabelsPath: csv file with the true labels (format: one column, no index)
  PredLabelsPath: csv file with the predicted labels (format: one column, no index)

  Returns
  -------
  Conf: confusion matrix
  MedF1 : median F1-score
  F1 : F1-score per class
  Acc : accuracy
  PercUnl : percentage of unlabeled cells
  PopSize : number of cells per cell type
  "
  # import python package: sklearn.metrics
  library(reticulate)
  # use_python('/home/drizzle_zhang/tools/anaconda3/bin/python3', required = T)
  use_python('/local/zy/tools/anaconda3/bin/python3', required = T)
  # py_config()
  py_module_available('sklearn')
  metrics <- import('sklearn.metrics')
  
  true_lab <- unlist(read.csv(TrueLabelsPath, stringsAsFactors = F))
  pred_lab <- unlist(read.csv(PredLabelsPath, stringsAsFactors = F))
  
  unique_true <- unlist(unique(true_lab))
  # unique_pred <- unlist(unique(pred_lab))
  
  # unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(true_lab,pred_lab)
  pop_size <- rowSums(conf)
  
  pred_lab = gsub('Node..','Node',pred_lab)
  
  conf_F1 <- table(true_lab,pred_lab)#,exclude = c('unassigned','Unassigned','Unknown','rand','Node','ambiguous'))
  
  F1 <- vector()
  sum_acc <- 0
  
  for (i in c(1:length(unique_true))) {
    findLabel = colnames(conf_F1) == (row.names(conf_F1))[i]
    if (sum(findLabel) > 0) {
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i,findLabel]
    } else {
      F1[i] = 0
    }
  }
  
  names(F1) <- names(pop_size)
  
  med_F1 <- median(F1)
  mean_F1 <- mean(F1)
  
  acc <- sum_acc/sum(conf_F1)
  
  total <- length(pred_lab)
  # num_unlab <- sum(pred_lab == 'unassigned') + sum(pred_lab == 'Unassigned') + sum(pred_lab == 'rand') + sum(pred_lab == 'Unknown') + sum(pred_lab == 'Node') + sum(pred_lab == 'ambiguous')
  num_unlab <- sum(!(pred_lab %in% unique_true))
  per_unlab <- num_unlab / total
  
  pred_lab[!(pred_lab %in% unique_true)] <- 'Unassigned'
  
  weighted.macro.F1 <- metrics$f1_score(true_lab, pred_lab, average = 'weighted')
  # metrics$f1_score(true_lab, pred_lab, average = 'macro')
  acc <- metrics$accuracy_score(true_lab, pred_lab)
  balanced_acc <- metrics$balanced_accuracy_score(true_lab, pred_lab)
  
  result <- list(Conf = conf, MedF1 = med_F1, F1 = F1, Mean_F1 = mean_F1, Acc = acc, 
                 balanced_acc = balanced_acc,
                 WMean_F1 = weighted.macro.F1,
                 PercUnl = per_unlab, PopSize = pop_size)
  
  return(result)
  
}
