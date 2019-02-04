# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


suppressPackageStartupMessages({
  library(scater)
  library(SingleCellExperiment)
  library(tidyverse)
  library(matrixStats)
  library(Rtsne)
  library(glue)
  library(scran)
  library(glmnet)
  library(plyr)
  library(magrittr)
  library(caret)
  library(doParallel)
  library(stats)
  library(pROC)
  library(ROCR)
  library(ggrepel)
  library(scales)
  
})








#gets genes with highest variance across cells in a singleCellExperiment object for the specified number of genes
variableGenes <- function(sce, marker_num){
  
  if(dim(sce)[1] > marker_num) {
    vars <- rowVars(as.matrix(logcounts(sce)))
    genes <- rownames(sce)
    df <- data.frame(genes, vars, stringsAsFactors = FALSE)
    df <- df[order(df$vars, decreasing=TRUE), ][1:marker_num,]
    
    return(sce[df$genes,])
    
    # topgenes <- arrange(df, desc(vars)) %>%
    #   head(n = numGenes) %>%
    #   .$genes
    #
    # return(sce[topgenes,])
  }
  else {
    return(sce)
  }
  
}








#Splits sce into testing and training data, but also filters the markers to reduce training time (if needed)
#Param: marker_num - number of markers to use
#Param: sce - single cell experiment with a clustering called mc_cluster
#Return: list(training_matrix, training_labels, testing_matrix, testing_labels)

sampleSplit <- function(sce, marker_num=2000) {
  
  if(dim(sce)[1] > marker_num) {
    sce<- variableGenes(sce, marker_num)
  }
  
  training_size <- ceiling(0.5 * length(colnames(sce)))
  
  training <- sce[, sample(colnames(sce), training_size)]
  testing <- sce[, setdiff(colnames(sce), colnames(training))]
  train_clus <- as.data.frame(as.character(training$clusters))
  test_clus <- as.data.frame(as.character(testing$clusters))
  
  
  training_matrix <- t(assay(training, "logcounts"))
  testing_matrix <- t(assay(testing, "logcounts"))
  training_label_df <- as.data.frame(training$clusters)
  training_labels <- as.data.frame(training$clusters)$`training$clusters`
  testing_labels <- as.data.frame(testing$clusters)$`testing$clusters`
  
  return(list(training_matrix, training_labels, testing_matrix, testing_labels))
  
}





cascer <- function(sce, clusters){
  
  sce$clusters <- clusters
  
  #split data
  data_splits <- sampleSplit(sce)
  train_y <- chartr("123456789", "ABCDEFGHI", data_splits[[2]])
  test_y <- chartr("123456789", "ABCDEFGHI", data_splits[[4]])
  
  #train glm
  trcntrl <- trainControl(method="cv", number=5, returnResamp="all",
               classProbs=TRUE,
               savePredictions = TRUE)

  fit <- train(data_splits[[1]], train_y, method = "glmnet",
                               trControl = trcntrl ,metric = "ROC",
                               tuneGrid = expand.grid(alpha = 1,
                                                      lambda = seq(0.001,0.1,by = 0.001)))

  
  
  probs <- predict(fit2, newdata = data_splits[[3]], type = "prob")
  
  classes <- chartr("ABCDEFGHI", "123456789", predict(fit2, newdata = data_splits[[3]], type = "raw")) %>%
    as.numeric()
  
  auc <- as.numeric(pROC::auc(multiclass.roc(as.numeric(classes), as.numeric(data_splits[[4]]))))
  
  obj <- list(predicted_classes=classes, auc=auc, response= probs, truths=data_splits[[4]])
  
  class(obj) <- "casc"
  
  return(obj)
}

multROC <- function(casc){
  truths <- casc$truths
  response <-casc$response
  #create list of rocs
  if(length(levels(truths) > 2)) {
    roc_l <- list()
    for(i in as.numeric(levels(as.factor(truths)))) {
      #binarize each level 
      single_t <- ifelse(truths == i, 1, 0)
      roc_l[i] <- list(pROC::roc(single_t, response[,i]))
      
    }
  }
  else {
    roc_l[1] <- list(pROC::roc(truths, response[,1]))
  }
  
  ggroc(roc_l)
}




