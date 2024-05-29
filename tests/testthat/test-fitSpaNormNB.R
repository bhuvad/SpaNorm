test_that("fitSpaNormNB internal function checks works", {
  set.seed(36)
  emat = matrix(rpois(6 * 4e3, 10), 6, 4e3)
  coords = cbind(1:4e3, 1:4e3)

  # maxit
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = 0), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = -1), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = 1.1), "integer")
})

test_that("fitNBGivenPsi internal function checks works", {
  Y = matrix(0, 6, 10)
  gmean = psi = rep(0, 6)
  W = matrix(0, 10, 3)
  alpha = matrix(0, 6, 3)

  # maxit
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = 0), "greater")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = -1), "greater")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = 1.1), "integer")

  # step.factor
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 0), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 1), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 1.1), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = -1.1), "0,1")
})

test_that("fitNBGivenPsi internal function checks works", {
  Y = matrix(0, 6, 10)
  gmean = psi = rep(0, 6)
  W = matrix(0, 10, 3)
  alpha = matrix(0, 6, 3)

  expect_error(checkNBParams(nrow(Y), ncol(Y), W[, -1], gmean, alpha, psi), "alpha.+W")
  expect_error(checkNBParams(nrow(Y), ncol(Y), W[-1, ], gmean, alpha, psi), "cells")
  expect_error(checkNBParams(nrow(Y), ncol(Y), W, gmean, alpha[, -1], psi), "alpha.+W")
  expect_error(checkNBParams(nrow(Y), ncol(Y), W, gmean, alpha[-1, ], psi), "genes")
  expect_error(checkNBParams(nrow(Y), ncol(Y), W, gmean[-1], alpha, psi), "genes")
  expect_error(checkNBParams(nrow(Y), ncol(Y), W, gmean[-1], alpha, psi[-1]), "genes")
  expect_equal(checkNBParams(nrow(Y), ncol(Y), W, gmean, alpha, psi), TRUE)
})