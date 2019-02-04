context("test-casc")

test_that("cascer returns correct types", {
  library(SingleCellExperiment)
  data(test_sce)
  data(test_casc_obj)

  casc_obj <- cascer(test_sce, test_sce$sc3_3_clusters)
  expect_is(casc_obj, "casc")

  plot <- multROC(test_casc_obj)
  expect_is(plot, "ggplot")

})
