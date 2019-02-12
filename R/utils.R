
#' @importFrom matrixStats rowVars
variableGenes <- function(sce, marker_num){
  
  if(dim(sce)[1] > marker_num) {
    vars <- rowVars(as.matrix(logcounts(sce)))
    genes <- rownames(sce)
    df <- data.frame(genes, vars, stringsAsFactors = FALSE)
    df <- df[order(df$vars, decreasing=TRUE), ][1:marker_num,]
    
    return(sce[df$genes,])
    
  }
  else {
    return(sce)
  }
  
}


sampleSplit <- function(sce, marker_num=2000) {
  
  if(dim(sce)[1] > marker_num) {
    sce<- variableGenes(sce, marker_num)
  }
  
  trainingSize <- ceiling(0.5 * length(colnames(sce)))
  
  training <- sce[, sample(colnames(sce), trainingSize)]
  testing <- sce[, setdiff(colnames(sce), colnames(training))]
  
  
  return(list(training,testing))
  
}


#' @importFrom caret trainControl train predict.train
#' @importFrom pROC auc multiclass.roc
#' @importFrom magrittr "%>%"
#' 
single_casc <- function(cluster_name, dataSplits){
  trainY <- chartr("123456789", "ABCDEFGHI", colData(dataSplits[[1]])[cluster_name][[1]])
  
  trcntrl <- trainControl(method="cv", number=5, returnResamp="all",
                          classProbs=TRUE,
                          savePredictions = TRUE)
  
  fit <- train(t(logcounts(dataSplits[[1]])), trainY, method = "glmnet",
               trControl=trcntrl, metric="Accuracy",
               tuneGrid=expand.grid(alpha=1,
                                    lambda=seq(0.001, 0.1, by=0.001)))
  
  
  probs <- predict.train(fit, newdata=t(logcounts(dataSplits[[2]])), type="prob")
  
  classes <- chartr("ABCDEFGHI", "123456789", predict(fit, newdata=t(logcounts(dataSplits[[2]])), type="raw")) %>%
    as.numeric()
  
  roc_l <- multROC(colData(dataSplits[[2]])[cluster_name][[1]], probs)
  auc <- avgAUC(roc_l)

  obj <- list(predicted_classes=classes, auc=auc, response= probs, truths=colData(dataSplits[[2]])[cluster_name][[1]] )
  
  class(obj) <- "casc"
  obj
}

#' Sample, Train, and Predict logistic regression model using singleCellExperiment and glmnet
#'
#' @param sce a singleCellExperiment
#' @param clusters A list of clusters, an array or list of integers of same length as number of cells in sce
#'
#' @return A casc object with predicted classes, auc, response, and truths
#' @import SingleCellExperiment
#' @export
cascer <- function(sce, clusters, marker_num=2000){
  

  k <- 0
  l <- list()
  for (c in clusters){
    k <- k + 1
    cluster_name <- paste0("casc_clusters_", k)
    colData(sce)[cluster_name] <- c
    l[k] <- cluster_name
  }

  dataSplits <- sampleSplit(sce, marker_num=marker_num)
  cascs <- lapply(l, single_casc, dataSplits=dataSplits)
  
  names(cascs) <- l
  
  class(cascs) <- "casc_list"
  
  
  cascs
  
}


#' Print a \code{cellassign} fit
#'
#' @param x An object of class \code{casc}
#' @param ... Additional arguments (unused)
#'
#'
#' @return Prints a structured representation of the \code{casc}
#' @import glue
#
#' @export
print.casc <- function(x, ...) {
  N <- length(levels(x$truths))
  
  cat(glue::glue("A casc object for {N} provided classes
                 ",
                 "Call x$predicted_classes to access predictions of class
                 ",
                 "Call x$auc to access area under the curve
                 ",
                 "Call x$response to access prediction response
                 ",
                 "Call x$truths to access original classes\n\n"))
}





#' calculate roc objects from a casc object, one for each cluster
#'
#' @param casc A casc object
#'
#' @return A list of roc objects
#' @importFrom  pROC roc
#' 
multROC <- function(truths, response){
  roc_l <- list()
  if(length(levels(truths)) > 2) {
    for(i in as.numeric(levels(as.factor(truths)))) {
      single_t <- ifelse(truths == i, 1, 0)
      roc_l[i] <- list(pROC::roc(single_t, response[,i]))
    }
  }
  else {
    roc_l[1] <- list(pROC::roc(truths, response[,1]))
  }

  roc_l
}

#' Plot multiple roc objects on the same ggplot
#'
#' @param casc A casc object
#'
#' @return A ggplot object with ROC curves plotted for each cluster
#' @importFrom  pROC ggroc
#' @importFrom ggplot2 theme element_blank
#' 
#' 
#' @export
multROCPlot <- function(casc){
  roc_l <- multROC(truths=casc$truths, response=casc$response)
  pROC::ggroc(roc_l) + theme(legend.title=element_blank())
}



#' Sample, Train, and Predict logistic regression model using singleCellExperiment and glmnet
#'
#' @param casc_list A casc_list object produced by cascer
#'
#' @return A ggplot scatterplot with AUCs plotted for each cluster
#' @importFrom ggplot2 ggplot geom_point theme_light
#' 
#' 
#' @export
aucPlot <- function(casc_list){
  
    df <- lapply(casc_list, function(x){x$auc}) %>%
      as.data.frame() %>%
      t()
    
    k <- lapply(casc_list, function(x){length(levels(x$truths))}) %>% 
      as.data.frame() %>%
      t()

    df %<>% cbind(k) %>%
      as.data.frame()
    
    colnames(df) <- c("AUC", "K")
    
    df$K <- as.factor(df$K)
    
    ggplot(df, aes(x=K, y=AUC)) + geom_point() + theme_light()
    
}


#' Sample, Train, and Predict logistic regression model using singleCellExperiment and glmnet
#'
#' @param roc_l A list of casc objects produced by cascer
#'
#' @return Average auc for the list, numeric
#' @importFrom roc auc
#' 
#' 
#' @export
avgAUC <- function(roc_l){
  auc_l <- lapply(roc_l, auc) %>% 
    lapply(as.numeric)
  mean(unlist(auc_l))
}





