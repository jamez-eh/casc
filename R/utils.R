#'
#'  Calculate most variable genes for SingleCellExperiment object
#'
#' @importFrom matrixStats rowVars
#'
#' @return sce with top marker_num most variable genes
#'
#' @keywords internal
variableGenes <- function(sce, marker_num) {
    if (dim(sce)[1] > marker_num) {
        vars <- rowVars(as.matrix(logcounts(sce)))
        genes <- rownames(sce)
        df <- data.frame(genes, vars, stringsAsFactors = FALSE)
        df <- df[order(df$vars, decreasing = TRUE),][seq_len(marker_num), ]
    
        return(sce[df$genes, ])
    }
    else {
        return(sce)
    }
}


#'
#'  Splits SingleCellExperiment object into test and training data
#'
#' @importFrom matrixStats rowVars
#'
#' @return a list, including both the training and test splits
#'
#' @keywords internal
sampleSplit <- function(sce, marker_num = 2000) {
    if (dim(sce)[1] > marker_num) {
        sce <- variableGenes(sce, marker_num)
    }
    trainingSize <- ceiling(0.5 * length(colnames(sce)))
    training <- sce[, sample(colnames(sce), trainingSize)]
    testing <- sce[, setdiff(colnames(sce), colnames(training))]
    
    return(list(training, testing))
}

#' calculate roc objects from a casc object, one for each cluster
#'
#' @param truths List of original class labels
#' @param response matrix or data.frame of probabilities of each class
#' @return A list of roc objects
#' @importFrom  pROC roc
multROC <- function(truths, response) {
    
    truths <- as.numeric(as.factor(truths))
    if (length(levels(as.factor(truths))) > 2) {
        roc_l <- lapply(seq_along(levels(as.factor(truths))), function(i){
            single_t <- ifelse(truths == i, 1, 0)
            pROC::roc(single_t, response[, i])
        })
    }
    else {
        roc_l <- list()
        roc_l[1] <- list(pROC::roc(truths, response[, 1]))
    }
    roc_l
}


#'
#'  Sample, Train, and Predict logistic regression model 
#'  using singleCellExperiment and glmnet, one object at a time
#'
#' @return a casc object
#'
#' @importFrom caret trainControl train predict.train
#' @importFrom magrittr "%>%"
#' @keywords internal
single_casc <- function(cluster_name, dataSplits, alpha = 0.5) {

    trainY <- colData(dataSplits[[1]])[[cluster_name]] %>%
        as.factor()
    
    original_levels <- levels(trainY)
    
    levels(trainY) <- do.call(paste0, 
                            replicate(7, sample(LETTERS, length(levels(trainY)), TRUE), 
                                                                                FALSE)) %>%
        sort()
    
    train_dat <- logcounts(dataSplits[[1]])[,!is.na(trainY)]
    trainY <- trainY[!is.na(trainY)]
    train_dat <- upsample(t(train_dat), trainY)
    
    trcntrl <- trainControl(
            method = "cv",
            number = 5,
            returnResamp = "all",
            classProbs = TRUE,
            savePredictions = TRUE
    )
    fit <-
        train(
        train_dat,
        trainY,
        method = "glmnet",
        trControl = trcntrl,
        metric = "Accuracy",
        tuneGrid = expand.grid(alpha = alpha,
                                lambda = seq(0.001, 0.1, by = 0.001))
    )
    
    testY <- colData(dataSplits[[2]])[[cluster_name]]
    test_dat <- logcounts(dataSplits[[2]])
    test_dat <- test_dat[,!is.na(testY)]
    testY <- testX[,!is.na(testY)]
    
    probs <-
        predict.train(fit, newdata = t(test_dat), 
                        type = "prob")

    classes <-predict(fit, newdata = t(test_dat), 
                type = "raw") %>%
        as.numeric() %>%
        as.factor()
    
    roc_l <- multROC(testY, probs) 
    auc <- avgAUC(roc_l)
    
    levels(classes) <- original_levels
    
    colnames(probs) <- original_levels
    
    obj <- list(
        predicted_classes = classes,
        auc = auc,
        response = probs,
        truths = as.factor(testY)
    )
    class(obj) <- "casc"
    obj
}


#' Sample, Train, and Predict logistic regression
#'  model using singleCellExperiment and glmnet.
#'
#' @param sce a singleCellExperiment
#' @param clusters A list of clusters, an array or list of 
#' integers of same length as number of cells in sce.
#' @param marker_num The top variable genes that 
#' will be used to filter the SCE by to reduce runtime.
#' @param alpha A parameter for logistic regression where 0 
#' is `ridge regression` and 1 is `lasso regression`.
#'
#' @return A list of casc objects with predicted classes, 
#' aucs, responses, and truths.
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' counts <- matrix(rnorm(40000, 10, 10), ncol=200, nrow=200)
#' sce <- SingleCellExperiment(assays = list(logcounts = counts))
#' colnames(sce) <- stringi::stri_rand_strings(200, 5)
#' rownames(sce) <- stringi::stri_rand_strings(200, 5)
#'
#' cluster_1 <- rep(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1), 20)
#' cluster_2 <- rep(c(0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 20)
#'
#' casc_list <- casc(sce, list(cluster_1, cluster_2), marker_num=1500)
#'
#'
#'
#' @import SingleCellExperiment
#' 
#' @importFrom magrittr "%>%"
#' @export
casc <- function(sce,
                clusters,
                marker_num = 2000,
                alpha = 0.5) {

    if(is.matrix(sce)){
        sce <- SingleCellExperiment(assays = list(logcounts = sce))
    }
    
    logcounts(sce) <- as.matrix(logcounts(sce))
    
    
    l <- lapply(seq_along(clusters), function(c) {
        cluster_name <- paste0("casc_", c)
        cluster_name
    })

    for(c in seq_along(l)){
        colData(sce)[l[[c]]] <- clusters[[c]]
    }
    
    dataSplits <- sampleSplit(sce, marker_num = marker_num)
    cascs <- lapply(l, single_casc, dataSplits = dataSplits, alpha = alpha)
    names(cascs) <- l
    class(cascs) <- "casc_list"
    
    cascs
}


#' Print a \code{casc}
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
    N <- dim(x$response)[2]
    cat(
        glue::glue(
        "A casc object for {N} provided classes
        ",
        "Call x$predicted_classes to access predictions of class
        ",
        "Call x$auc to access area under the curve
        ",
        "Call x$response to access prediction response
        ",
        "Call x$truths to access original classes\n\n"
        )
    )
}


#' Plot multiple roc objects on the same ggplot
#'
#' @param casc A casc object
#'
#' @return A ggplot object with ROC curves plotted for each cluster
#' @importFrom  pROC ggroc
#' @importFrom ggplot2 theme element_blank
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' counts <- matrix(rnorm(40000, 10, 10), ncol=200, nrow=200)
#' sce <- SingleCellExperiment(assays = list(logcounts = counts))
#' colnames(sce) <- stringi::stri_rand_strings(200, 5)
#' rownames(sce) <- stringi::stri_rand_strings(200, 5)
#'
#' cluster_1 <- rep(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1), 20)
#' cluster_2 <- rep(c(0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 20)
#'
#' casc_list <- casc(sce, list(cluster_1, cluster_2), marker_num=1500)
#' multROCPlot(casc_list$casc_1)
#'
#' @export
multROCPlot <- function(casc) {
    roc_l <- multROC(truths = casc$truths, response = casc$response)
    pROC::ggroc(roc_l) + theme(legend.title = element_blank())
}


#' Sample, Train, and Predict logistic regression 
#' model using singleCellExperiment and glmnet
#'
#' @param casc_list A casc_list object produced by casc
#'
#' @return A ggplot scatterplot with AUCs plotted for each cluster
#' @importFrom ggplot2 ggplot geom_point theme_light
#' @importFrom magrittr "%<>%"
#' @importFrom ggplot2 aes
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' counts <- matrix(rnorm(40000, 10, 10), ncol=200, nrow=200)
#' sce <- SingleCellExperiment(assays = list(logcounts = counts))
#' colnames(sce) <- stringi::stri_rand_strings(200, 5)
#' rownames(sce) <- stringi::stri_rand_strings(200, 5)
#'
#' cluster_1 <- rep(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1), 20)
#' cluster_2 <- rep(c(0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 20)
#'
#' casc_list <- casc(sce, list(cluster_1, cluster_2), marker_num=1500)
#' aucPlot(casc_list)
#'
#' @export
aucPlot <- function(casc_list) {
    df <- t(as.data.frame(lapply(casc_list, function(x) {
        x$auc
    })))

    k <- t(as.data.frame(lapply(casc_list, function(x) {
        length(x[[3]])
    })))
    
    df <- as.data.frame(cbind(df, k)) 
    colnames(df) <- c("AUC", "K")
    df$K <- as.factor(df$K)
    
    ggplot(df, aes(x = K, y = AUC)) + geom_point() + theme_light()
}


#' Finds the average AUC from a list of roc objects
#'
#' @param roc_l A list of casc objects produced by casc
#' @importFrom pROC auc
#' @importFrom magrittr "%>%"
#' @return Average auc for the list, numeric
#'
avgAUC <- function(roc_l) {
    auc_l <- lapply(roc_l, auc) %>%
        lapply(as.numeric)
    mean(unlist(auc_l))
}