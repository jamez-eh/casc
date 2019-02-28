context("test-casc")

test_that("cascer returns correct types", {
  library(SingleCellExperiment)
  data(test_sce)
  data(test_casc_list)

  casc_obj <- cascer(test_sce, list(test_sce$sc3_3_clusters, test_sce$sc3_2_clusters))
  
  expect_equal(length(casc_obj), 2)
  
  expect_is(casc_obj$casc_clusters_1, "casc")
  
  roc_l <- multROC(casc_obj$casc_clusters_1$truths, casc_obj$casc_clusters_1$response)
  expect_is(roc_l[[1]], "roc")
  
  
  plot <- multROCPlot(test_casc_list$casc_clusters_1)
  expect_is(plot, "ggplot")

  plot_2 <- aucPlot(roc_l)
  expect_is(plot_2, "ggplot")
  
  
})

