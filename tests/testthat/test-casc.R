context("test-casc")

test_that("cascer returns correct types", {
  library(SingleCellExperiment)
  data(test_sce)
  data(test_casc_list)

  casc_obj <- cascer(test_sce, test_sce$sc3_3_clusters, test_sce$sc3_2_clusters)
  
  expect_equal(length(casc_obj), 2)
  
  expect_is(casc_obj$casc_clusters_1, "casc")
  

  plot <- multROC(test_casc_list$casc_clusters_1)
  expect_is(plot, "ggplot")

})

