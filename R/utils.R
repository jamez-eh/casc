

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

single_casc <- function(cluster_name, dataSplits){
  trainY <- chartr("123456789", "ABCDEFGHI", colData(dataSplits[[1]])[cluster_name][[1]])
  
  trcntrl <- trainControl(method="cv", number=5, returnResamp="all",
                          classProbs=TRUE,
                          savePredictions = TRUE)
  
  fit <- train(t(logcounts(dataSplits[[1]])), trainY, method = "glmnet",
               trControl=trcntrl, metric="Accuracy",
               tuneGrid=expand.grid(alpha=1,
                                    lambda=seq(0.001, 0.1, by=0.001)))
  
  
  probs <- predict(fit, newdata=t(logcounts(dataSplits[[2]])), type="prob")
  
  classes <- chartr("ABCDEFGHI", "123456789", predict(fit, newdata=t(logcounts(dataSplits[[2]])), type="raw")) %>%
    as.numeric()
  
  auc <- as.numeric(pROC::auc(multiclass.roc(as.numeric(classes), as.numeric(colData(dataSplits[[2]])[cluster_name][[1]]))))
  
  obj <- list(predicted_classes=classes, auc=auc, response= probs, truths=colData(dataSplits[[2]])[cluster_name][[1]] )
  
  class(obj) <- "casc"
  obj
}

#' Sample, Train, and Predict logistic regression model using singleCellExperiment and glmnet
#'
#' @param sce a singleCellExperiment
#' @param clusters A cluster, an array or list of integers of same length as number of cells in sce
#'
#' @return A casc object with predicted classes, auc, response, and truths
#'
#' @export
cascer <- function(sce, ...){
  
  clusters <- list(...)
  
  k <- 0
  l <- list()
  for (c in clusters){
    k <- k + 1
    cluster_name <- paste0("casc_clusters_", k)
    colData(sce)[cluster_name] <- c
    l[k] <- cluster_name
  }
  sce$clusters <- clusters
  
  dataSplits <- sampleSplit(sce)
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





#' Sample, Train, and Predict logistic regression model using singleCellExperiment and glmnet
#'
#' @param casc A casc object produced by cascer
#'
#' @return A ggplot object with ROC curves plotted for each cluster
#'
#' @export
multROC <- function(casc){
  truths <- casc$truths
  response <-casc$response
  if(length(levels(truths) > 2)) {
    roc_l <- list()
    for(i in as.numeric(levels(as.factor(truths)))) {
      single_t <- ifelse(truths == i, 1, 0)
      roc_l[i] <- list(pROC::roc(single_t, response[,i]))
    }
  }
  else {
    roc_l[1] <- list(pROC::roc(truths, response[,1]))
  }

  pROC::ggroc(roc_l) + theme(legend.title=element_blank())
}




