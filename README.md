#  casc

`casc` trains a logisitic regression model on provided single-cell RNAseq clusters. It creates an ROC curve for each cluster vs the others to help decide the appropriate number of clusters.

## Installation

`casc` can be installed via `devtools` and github:

``` r
install.packages("devtools") 
devtools::install_github("jamez-eh/casc")
```
## Useage

`casc` has the following parameters:

* `sce`: A `SingleCellExperiment` with normalized logcounts as an assay or a logcounts matrix with cells as columns and genes as rows.
* `clusters`: A list of clusterings to evaluate.
* `alpha`: A parameter for logistic regression. alpha = 1 represents the lasso penalty and alpha = 0 represents the ridge penalty.

#Author

James Hopkins