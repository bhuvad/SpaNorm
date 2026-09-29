# tpsBasis(): the exported spatial basis. On its own coordinates it must be
# bs.tps() exactly -- the basis SpaNorm's fit is built from -- and built on a
# reference it must be the same function of position wherever it is evaluated.

test_that("tpsBasis on its own coordinates is bs.tps(), bit for bit", {
  set.seed(1)
  # a near-square field, a 3:1 rectangle (the aspect-ratio rule gives the
  # short axis fewer degrees of freedom) and SpaNorm's own unit-range scaling
  x <- runif(300); y <- runif(300)
  xr <- c(0, 3, runif(248, 0, 3)); yr <- c(0, 1, runif(248, 0, 1))
  sx <- (xr - min(xr)) / (max(xr) - min(xr)) - 0.5
  sy <- (yr - min(yr)) / (max(yr) - min(yr)) - 0.5
  for (df in c(1L, 2L, 3L, 6L)) {
    expect_identical(tpsBasis(x, y, df), bs.tps(x, y, df))
    expect_identical(tpsBasis(xr, yr, df), bs.tps(xr, yr, df))
    expect_identical(tpsBasis(sx, sy, df), bs.tps(sx, sy, df))
  }
  expect_equal(attr(tpsBasis(xr, yr, 6L), "df.tps"), c(6, 2))
})

test_that("tpsBasis reproduces the biology basis of a real SpaNorm fit", {
  spe <- .polish_spe()
  fit <- S4Vectors::metadata(spe)$SpaNorm
  coords <- SpatialExperiment::spatialCoords(spe)
  # fitSpaNorm() scales each axis to [-0.5, 0.5] before building the basis
  sc <- apply(coords, 2, function(v) (v - min(v)) / (max(v) - min(v)) - 0.5)
  B <- tpsBasis(sc[, 1], sc[, 2], df = max(fit$df.tps[1:2]))
  expect_equal(unname(fit$W[, fit$wtype == "biology"]), unname(B[, ]),
               tolerance = 0, ignore_attr = TRUE)
})

test_that("built on a reference, the basis is one function of position", {
  set.seed(2)
  x <- runif(400, 0, 2); y <- runif(400)
  full <- tpsBasis(x, y, df = 4)
  sub <- which(x < 0.7 & y > 0.2)
  part <- tpsBasis(x[sub], y[sub], df = 4, ref.x = x, ref.y = y)
  # the same knots and the same centring, so the subset's rows are the
  # reference basis's rows
  expect_equal(unclass(part)[, ], unclass(full)[sub, ], tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_identical(attr(part, "scaled:center"), attr(full, "scaled:center"))
  expect_identical(attr(part, "df.tps"), attr(full, "df.tps"))
  # centred over the REFERENCE: zero column means there, not on the subset
  expect_equal(unname(colMeans(full)), rep(0, ncol(full)), tolerance = 1e-12)
  expect_gt(max(abs(colMeans(part))), 1e-3)
  # a point outside the reference range is extrapolated, not refused
  out <- tpsBasis(c(-0.5, 2.5), c(0.5, 1.2), df = 4, ref.x = x, ref.y = y)
  expect_true(all(is.finite(out)))
})

test_that("df can be given per axis, and centring can be turned off", {
  set.seed(3)
  x <- runif(200, 0, 5); y <- runif(200)
  B <- tpsBasis(x, y, df = c(3, 3))
  expect_equal(ncol(B), 9L)
  expect_equal(attr(B, "df.tps"), c(3L, 3L))
  U <- tpsBasis(x, y, df = c(3, 3), center = FALSE)
  expect_null(attr(U, "scaled:center"))
  expect_equal(sweep(unclass(U)[, ], 2L, colMeans(U)), unclass(B)[, ],
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("tpsBasis validates its inputs", {
  expect_error(tpsBasis(1:3, 1:4), "same length")
  expect_error(tpsBasis(c(1, NA, 3), 1:3), "finite")
  expect_error(tpsBasis(1:5, 1:5, df = 0), "positive integers")
  expect_error(tpsBasis(1:5, 1:5, df = 2.5), "positive integers")
  expect_error(tpsBasis(1:5, 1:5, df = c(1, 2, 3)), "positive integers")
  expect_error(tpsBasis(1:5, 1:5, ref.x = 1:3, ref.y = 1:4), "same length")
  expect_error(tpsBasis(1:5, rep(1, 5), df = 3), "no range")
})
