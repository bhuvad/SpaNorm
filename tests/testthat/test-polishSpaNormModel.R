mk_fit <- function(shared = TRUE) {
  set.seed(5)
  n <- 50; G <- 4
  W <- cbind(logLS = rnorm(n), b1 = rnorm(n), b2 = rnorm(n), l2 = rnorm(n))
  al <- matrix(rnorm(G * 4), G, 4)
  al[, 1] <- if (shared) 1.02 else seq(0.9, 1.1, length.out = G)
  SpaNormFit(ngenes = G, ncells = n, gene.model = "nb", df.tps = 2L,
             sample.p = 0.5, lambda.a = c(1e-4, 2e-4), batch = NULL, W = W,
             alpha = al, gmean = rnorm(G), psi = rep(0.2, G),
             wtype = factor(c("ls", "biology", "biology", "ls")),
             loglik = 0, sampling = factor(rep(c("glm", "all"), each = n / 2)))
}

test_that("the problem maps SpaNorm's model: intercept, offset, wtype penalty", {
  f <- mk_fit()
  p <- .spaNormPolishProblem(f)
  expect_equal(unname(p$X[, 1]), rep(1, 50))
  expect_equal(p$pen, c(0, 1e-4 * 50, 1e-4 * 50, 2e-4 * 50))
  expect_equal(p$offset, 1.02 * f$W[, 1])
  expect_equal(p$A0, cbind(f$gmean, f$alpha[, -1]))
  expect_equal(sum(.spaNormPolishProblem(f, "fit")$cells_idx), 25)
})

test_that("a per-gene library-size coefficient is refused", {
  expect_error(.spaNormPolishProblem(mk_fit(shared = FALSE)), "not .*shared")
})

test_that("column 1 of X is named '(gmean)' even when W has no column names", {
  f <- mk_fit()
  f@W <- unname(f@W)
  expect_null(colnames(f@W))
  p <- .spaNormPolishProblem(f)
  expect_identical(colnames(p$X)[1], "(gmean)")
})

test_that(".spaNormPolishProblem(cells = 'fit') selects glm/dispersion cells, not 'all'", {
  f <- mk_fit()
  p <- .spaNormPolishProblem(f, "fit")
  expect_identical(p$cells_idx, f$sampling != "all")
  expect_identical(nrow(p$X), sum(f$sampling != "all"))
})

test_that("a fit saved before the polish slot existed reads as unpolished", {
  f <- mk_fit()
  expect_false(isPolished(f))
  # simulate an old object: drop the slot from the attributes. Ruling 3 notes
  # this trick alone is not a faithful stand-in for a real pre-slot RDS (see
  # the fixture-based test below), so it is kept as a cheap extra check, not
  # the only one.
  old <- f; attr(old, "polish") <- NULL
  expect_false(isPolished(old))
  expect_identical(.polishSlot(old), list())
})

test_that("an RDS written before the polish slot existed reads as unpolished", {
  # Generated once, under the INSTALLED pre-slot SpaNorm 1.7.13 (the dev
  # library, before this branch added the `polish` slot) -- not
  # devtools::load_all(), which would pick up the in-progress working tree:
  #
  #   R_LIBS=/scratch/user/uqdbhuva/Rlib_polishmove Rscript -e '
  #     library(SpaNorm)
  #     stopifnot(packageVersion("SpaNorm") == "1.7.13")
  #     stopifnot(!"polish" %in% methods::slotNames("SpaNormFit"))
  #     mk_fit <- function(shared = TRUE) {
  #       set.seed(5)
  #       n <- 50; G <- 4
  #       W <- cbind(logLS = rnorm(n), b1 = rnorm(n), b2 = rnorm(n), l2 = rnorm(n))
  #       al <- matrix(rnorm(G * 4), G, 4)
  #       al[, 1] <- if (shared) 1.02 else seq(0.9, 1.1, length.out = G)
  #       SpaNorm:::SpaNormFit(ngenes = G, ncells = n, gene.model = "nb", df.tps = 2L,
  #         sample.p = 0.5, lambda.a = c(1e-4, 2e-4), batch = NULL, W = W,
  #         alpha = al, gmean = rnorm(G), psi = rep(0.2, G),
  #         wtype = factor(c("ls", "biology", "biology", "ls")),
  #         loglik = 0, sampling = factor(rep(c("glm", "all"), each = n / 2)))
  #     }
  #     saveRDS(mk_fit(), "spanormfit_1713_noslot.rds", version = 2)'
  #
  # (SpaNormFit() is not exported by SpaNorm, so the generating script reaches
  # it via SpaNorm:::SpaNormFit().)
  old <- readRDS(testthat::test_path("fixtures", "spanormfit_1713_noslot.rds"))
  expect_false(methods::.hasSlot(old, "polish"))
  expect_false(isPolished(old))
  expect_identical(.polishSlot(old), list())
})

test_that("validObject() on the pre-slot fixture errors (documented R S4 behaviour)", {
  # R's S4 machinery checks slot presence against the CURRENT class
  # definition before validSpaNormFit() (registered via setValidity()) ever
  # runs, so an object serialized under the pre-slot class does not validate
  # under the new one -- observed message:
  #   invalid class "SpaNormFit" object: slots in class definition but not
  #   in object: "polish"
  # `old@polish` errors the same way ("no slot of name 'polish' for this
  # object"). isPolished()/.polishSlot() (and `$`, which is
  # slot(x, name)) never call validObject() or read @polish directly, so
  # none of this blocks the "reads as unpolished" behaviour tested above --
  # it is a property of validObject()/@ alone.
  old <- readRDS(testthat::test_path("fixtures", "spanormfit_1713_noslot.rds"))
  expect_error(validObject(old), "slots in class definition but not in object")
  expect_error(old@polish, "no slot of name")
})

test_that("validity: polish must be an empty list, or have 'settings' and 'genes'", {
  f <- mk_fit()
  expect_true(isVirtualClass("SpaNormFit") || is(f, "SpaNormFit"))

  f@polish <- list(settings = list(psi.method = "fixed"), genes = data.frame())
  expect_true(validObject(f))

  f@polish <- list(settings = list(psi.method = "fixed")) # no 'genes'
  expect_error(validObject(f), "polish")

  f@polish <- list(nonsense = 1) # neither element
  expect_error(validObject(f), "polish")
})

test_that("show() adds one line when polished, and nothing extra when not", {
  f <- mk_fit()
  out_unpolished <- capture.output(show(f))
  expect_false(any(grepl("^polished:", out_unpolished)))

  f@polish <- list(
    settings = list(psi.method = "fixed", ls = "fixed", cells = "all"),
    genes = data.frame(gene = character(0))
  )
  out_polished <- capture.output(show(f))
  expect_true(any(out_polished ==
    "polished: psi.method = fixed, ls = fixed, cells = all"))
  # everything else about the display is unchanged
  expect_identical(setdiff(out_unpolished, out_polished), character(0))
})
