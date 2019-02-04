context("test-casc")

test_that("cascer returns correct types", {
  test_sce <- data(test_sce)


  casc_obj <- cascer(test_sce, sce$sc3_3_clusters)
  expect_is(casc_obj, "casc")
})
