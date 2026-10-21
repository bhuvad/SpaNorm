test_that("aesColnames() returns names without backticks", {
  f = function(...) aesColnames(rlang::enquos(...))
  expect_identical(f(colour = !!rlang::sym("HLA-B"), fill = EPCAM), c("HLA-B", "EPCAM"))
  expect_identical(f(), character(0))
})

test_that("plotSpatial() plots genes with non-syntactic names", {
  data("HumanDLPFC", package = "SpaNorm", envir = environment())
  spe = HumanDLPFC
  rownames(spe)[1] = "HLA-B"
  p = plotSpatial(spe, what = "expression", assay = "counts",
                  colour = !!rlang::sym("HLA-B"), size = 0.2)
  expect_s3_class(p, "ggplot")
  expect_true("HLA-B" %in% colnames(p$layers[[1]]$data))
  expect_no_error(ggplot2::ggplot_build(p))
})
