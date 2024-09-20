test_that("checks of model parameters work", {
  ngenes = 10
  ncells = 5
  nW = 3

  params = list(
    ngenes,
    ncells,
    nW,
    gene.model = "nb",
    df.tps = 6L,
    sample.p = 1,
    lambda.a = 0.5,
    batch = NULL,
    W = matrix(0, ncells, nW),
    alpha = matrix(0, ngenes, nW),
    gmean = rep(0, ngenes),
    psi = rep(0, ngenes),
    isbio = rep(TRUE, nW),
    loglik = rep(0, 5)
  )

  expect_equal(is(do.call(SpaNormFit, params))[[1]], "SpaNormFit")

  # df.tps
  tmp = params; tmp$df.tps = 0L
  expect_error(do.call(SpaNormFit, tmp), "greater")

  # gene.model
  tmp = params; tmp$gene.model = "gamma"
  expect_error(do.call(SpaNormFit, tmp), "nb")

  # sample.p
  tmp = params; tmp$sample.p = 0
  expect_error(do.call(SpaNormFit, tmp), "interval")
  tmp = params; tmp$sample.p = 1.1
  expect_error(do.call(SpaNormFit, tmp), "interval")
  tmp = params; tmp$sample.p = -1
  expect_error(do.call(SpaNormFit, tmp), "interval")

  # lambda.a
  tmp = params; tmp$lambda.a = -1
  expect_error(do.call(SpaNormFit, tmp), "greater")

  # loglik
  tmp = params; tmp$loglik[1] = 1
  expect_error(do.call(SpaNormFit, tmp), "less")
})

test_that("checks of model values work", {
  ngenes = 10
  ncells = 5
  nW = 3
  nbatch = 2

  params = list(
    ngenes,
    ncells,
    nW,
    gene.model = "nb",
    df.tps = 6L,
    sample.p = 1,
    lambda.a = 0.5,
    batch = NULL,
    W = matrix(0, ncells, nW),
    alpha = matrix(0, ngenes, nW),
    gmean = rep(0, ngenes),
    psi = rep(0, ngenes),
    isbio = rep(TRUE, nW),
    loglik = rep(0, 5)
  )

  expect_equal(is(do.call(SpaNormFit, params))[[1]], "SpaNormFit")

  # W
  tmp = params; tmp$W = tmp$W[, -1]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$W = tmp$W[-1, ]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$W = tmp$W[1, ]
  expect_error(do.call(SpaNormFit, tmp), "class")
  tmp = params; tmp$W = tmp$W[, 1]
  expect_error(do.call(SpaNormFit, tmp), "class")
  tmp = params; tmp$W[1, 1] = NA
  expect_error(do.call(SpaNormFit, tmp), "missing")

  # alpha
  tmp = params; tmp$alpha = tmp$alpha[, -1]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$alpha = tmp$alpha[-1, ]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$alpha = tmp$alpha[1, ]
  expect_error(do.call(SpaNormFit, tmp), "class")
  tmp = params; tmp$alpha = tmp$alpha[, 1]
  expect_error(do.call(SpaNormFit, tmp), "class")
  tmp = params; tmp$alpha[1, 1] = NA
  expect_error(do.call(SpaNormFit, tmp), "missing")

  # gmean
  tmp = params; tmp$gmean = tmp$gmean[-1]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$gmean = tmp$gmean[1]
  expect_error(do.call(SpaNormFit, tmp), "match")
  tmp = params; tmp$gmean[1] = NA
  expect_error(do.call(SpaNormFit, tmp), "missing")

  # psi
  tmp = params; tmp$psi = tmp$psi[-1]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$psi = tmp$psi[1]
  expect_error(do.call(SpaNormFit, tmp), "match")
  tmp = params; tmp$psi[1] = NA
  expect_error(do.call(SpaNormFit, tmp), "missing")

  # isbio
  tmp = params; tmp$isbio = tmp$isbio[-1]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$isbio = tmp$isbio[1]
  expect_error(do.call(SpaNormFit, tmp), "match")
  tmp = params; tmp$isbio[1] = NA
  expect_error(do.call(SpaNormFit, tmp), "missing")
  tmp = params; tmp$isbio = !tmp$isbio
  expect_error(do.call(SpaNormFit, tmp), "at least")

  # batch - vector
  params$batch = rep("A", ncells)
  expect_equal(is(do.call(SpaNormFit, params))[[1]], "SpaNormFit")
  tmp = params; tmp$batch = tmp$batch[-1]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$batch = tmp$batch[1]
  expect_error(do.call(SpaNormFit, tmp), "match")
  tmp = params; tmp$batch[1] = NA
  expect_error(do.call(SpaNormFit, tmp), "missing")
  
  # batch - matrix
  params$batch = matrix("A", ncells, nbatch)
  expect_equal(is(do.call(SpaNormFit, params))[[1]], "SpaNormFit")
  tmp = params; tmp$batch = tmp$batch[-1, ]
  expect_error(do.call(SpaNormFit, tmp), "not match")
  tmp = params; tmp$batch[1, 1] = NA
  expect_error(do.call(SpaNormFit, tmp), "missing")
})

