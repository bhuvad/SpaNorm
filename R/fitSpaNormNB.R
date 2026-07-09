#' @importFrom stats mad median quantile model.matrix pnbinom dnbinom qnbinom
NULL

fitSpaNormNB <- function(Y, W, idx, maxit.psi = 25, tol = 1e-4, maxn.psi = 500, ..., backend = c("auto", "cpu", "gpu"), msgfun = message) {
  backend = match.arg(backend)
  # parameter checks
  if (maxit.psi <= 0) {
    stop("'maxit.psi' should be greater than 0")
  }
  if (maxit.psi - as.integer(maxit.psi) != 0) {
    stop("'maxit.psi' should be integers")
  }

  # message function for inner loop
  msgfunNB = \(x, ...) {
    \(...) msgfun(x, ...)
  }

  # subset data points used for model fitting
  Ysub = toGPUMatrix(as.matrix(Y[, idx, drop = FALSE]), backend = backend)
  Wsub = toGPUMatrix(W[idx, , drop = FALSE], backend = backend)
  nW = ncol(Wsub)
  nsub = sum(idx)

  # setting initial regression parameter estimates for all genes
  gmean = rowMeans_gpu(log(Ysub + 1)) # rowMeans(Z)
  alpha = matrix(0, nrow(Ysub), ncol(W))
  # alpha for logLS
  alpha[, 1] = 1
  alpha = toGPUMatrix(alpha, backend = backend)

  # subset to a maximum of 50 cells/spots for dispersion parameter estimation (for speed-up)
  nsub.psi = min(maxn.psi, nsub)
  psi.idx = rep(FALSE, nsub)
  psi.idx[sample.int(nsub, size = nsub.psi)] = TRUE
  psi.cols = which(psi.idx)
  # realise only the dispersion subset (genes x nsub.psi), not all sampled cells
  Ysub.psi = as.matrix(if (is_torch_tensor(Ysub)) Ysub[, psi.cols] else Ysub[, psi.cols, drop = FALSE])
  psi = rep(0, nrow(Ysub))
  
  # convergence trackers
  conv = FALSE
  logl.psi = c()
  iter = 1

  # map sample ids used for estimation for export
  idx = as.numeric(idx)
  idx[idx == 1][psi.idx] = 2
  idx = as.factor(c("all", "glm", "dispersion")[idx + 1])

  # outer loop of IRLS for estimating disp parameter
  while (!conv) {
    msgfun(sprintf("iter:%3d, estimating gene-wise dispersion", iter))
    # calculate/update dispersion parameter (psi) estimates for all genes
    Wa = tcrossprod_gpu(alpha, Wsub) # offsets
    # slice the dispersion columns before realising (avoids materialising all
    # sampled cells of Wa on the CPU each outer iteration)
    offs.psi = as.matrix(if (is_torch_tensor(Wa)) Wa[, psi.cols] else Wa[, psi.cols, drop = FALSE])
    psi.tmp = tryCatch(
      {
        edgeR::estimateDisp(Ysub.psi, as.matrix(rep(1, nsub.psi)), offset = offs.psi, tagwise = TRUE, robust = TRUE)$tagwise.dispersion
      },
      error = function(e) {
        rep(NA, length(psi))
      }
    )
    valid.psi = !is.na(psi.tmp) & !is.infinite(psi.tmp)
    psi[valid.psi] = psi.tmp[valid.psi]

    # calculate initial loglik for this iteration
    loglik = as.numeric(colSums_gpu(dnbinom_gpu(Ysub, mu = exp(add_vec_mat_gpu(gmean, Wa, backend = backend)), size = 1/psi, log = TRUE)))
    msgfun(sprintf("iter:%3d, log-likelihood: %f", iter, sum(loglik)))

    # fit NB given dispersion estimates (psi) and extract required components
    msgfun(sprintf("iter:%3d, fitting NB model", iter))
    fit.nb = fitNBGivenPsi(Ysub, Wsub, psi, ..., gmean = gmean, alpha = alpha, loglik = loglik, backend = backend, msgfun = msgfunNB(sprintf("iter:%3d, ", iter)))
    gmean = fit.nb$gmean
    alpha = fit.nb$alpha
    loglik = fit.nb$loglik

    # check convergence
    logl.psi = c(sum(loglik), logl.psi) # push to stack
    conv.logl = FALSE
    iter = iter + 1 # similar post-increment to fitNBGivenPsi
    if (iter > 2) { # should be 2 because iter is post-incremented
      conv.logl = ifelse((logl.psi[1] - logl.psi[2]) / abs(logl.psi[2]) < tol, TRUE, FALSE)
    }
    conv = conv.logl | iter > maxit.psi
  }

  # convergence message
  if (iter > maxit.psi) {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (did not converged)", iter, logl.psi[1]))
  } else {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (converged)", iter, logl.psi[1]))
  }

  # final values
  fit.spanorm = list(
    gmean = as.vector(gmean),
    alpha = as.matrix(alpha),
    psi = as.vector(psi),
    sampling = idx,
    loglik = rev(logl.psi)
  )

  return(fit.spanorm)
}

fitNBGivenPsi <- function(Ysub, Wsub, psi, lambda.a, gmean = NULL, alpha = NULL, step.factor = 0.5, maxit.nb = 50, tol = 1e-4, loglik = NULL, backend = c("auto", "cpu", "gpu"), is.spanorm = FALSE, msgfun = message) {
  backend = match.arg(backend)
  # parameter checks
  if (any(lambda.a < 0)) {
    stop("'lambda.a' should be positive")
  }
  if (!length(lambda.a) %in% c(1, ncol(Wsub) - 1)) {
    stop("'lambda.a' should either be a single value or a vector of length equal to the number of columns in W - 1")
  }
  if (length(lambda.a) == 1) {
    lambda.a = rep(lambda.a, ncol(Wsub) - 1)
  }
  if (step.factor <= 0 | step.factor >= 1) {
    stop("'step.factor' should be in the interval (0,1)")
  }
  if (maxit.nb <= 0) {
    stop("'maxit.nb' should be greater than 0")
  }
  if (maxit.nb - as.integer(maxit.nb) != 0) {
    stop("'maxit.nb' should be integers")
  }

  # convert matrices to GPU matrices if available
  Ysub = toGPUMatrix(Ysub, backend = backend)
  Wsub = toGPUMatrix(Wsub, backend = backend)
  lambda.a = diag_mat(lambda.a, backend = backend)
  psi = toGPUVector(psi, backend = backend)

  if (is.null(gmean)) {
    gmean = rowMeans_gpu(log(Ysub + 1)) # rowMeans(Z)
  }
  gmean = toGPUVector(gmean, backend = backend)
  if (is.null(alpha)) {
    alpha = matrix(0, nrow(Ysub), ncol(Wsub))
    # alpha for logLS
    alpha[, 1] = 1
  }
  alpha = toGPUMatrix(alpha, backend = backend)

  # set the IRLS step size to 1
  nsub = ncol(Ysub)
  nW = ncol(Wsub)
  step = rep(1, nsub)
  
  if (is.null(loglik))  {
    # if not pre-computed (e.g., by outer loop), compute
    # helpful pattern for external use of the fitting function
    lmu.hat = add_vec_mat_gpu(gmean, tcrossprod_gpu(alpha, Wsub), backend = backend)
    # psi (length n_genes == nrow) broadcasts per-row, matching the outer()
    # form used previously but without materialising an n_genes x nsub matrix
    # (outer() also fails on a torch-tensor psi)
    sig.inv = 1 / add_vec_mat_gpu(psi, exp(-lmu.hat), backend = backend)
    loglik = as.numeric(colSums_gpu(dnbinom_gpu(Ysub, mu = exp(lmu.hat), size = 1 / psi, log = TRUE)))
  }

  # convergence trackers
  conv = FALSE
  logl.beta = c(sum(loglik)) # stack of loglik values - push values
  iter = 1

  # initialise halving counters
  halving = 0

  while (!conv) {
    # message
    msgfun(sprintf("iter:%3d, log-likelihood: %f", iter, logl.beta[1]))

    # saving new best estimate (in case we need to go back in our search)
    best.gmean = copy(gmean)
    best.a = copy(alpha)

    # calc working vector Z for a subset of cells
    lmu.hat = add_vec_mat_gpu(gmean, tcrossprod_gpu(alpha, Wsub), backend = backend)
    # working residual, scaled per-cell (per-column) by the IRLS step vector.
    # sweep() cannot run on a torch tensor (it mixes an R vector, which lands on
    # the CPU, with a device tensor), so scale explicitly on the tensor path
    # while keeping the exact base-R behaviour on the CPU path.
    Zr = (Ysub + 0.01) / (exp(lmu.hat) + 0.01) - 1
    if (is_torch_tensor(Zr)) {
      step.t = torch::torch_tensor(as.numeric(step), dtype = Zr$dtype, device = Zr$device)$unsqueeze(1L)
      Z = lmu.hat + Zr * step.t
    } else {
      Z = lmu.hat + sweep(Zr, 2, step, "*")
    }

    # store current alpha est
    alpha.old = copy(alpha)

    # update alpha for all genes
    sig.inv = 1 / add_vec_mat_gpu(psi, exp(-lmu.hat), backend = backend)
    wt.cell = colMeans_gpu(sig.inv)
    # prevent outlier large weight; stay on-device for tensors to avoid a
    # per-iteration GPU->CPU->GPU round trip through quantile()
    if (is_torch_tensor(wt.cell)) {
      wt.cell = torch::torch_minimum(wt.cell, torch::torch_quantile(wt.cell, 0.98))
    } else {
      wt.cell = pmin(wt.cell, quantile(wt.cell, probs = 0.98))
    }
    Wsub.wt = mult_vec_mat_gpu(wt.cell, Wsub) # row-scale W by cell weights (avoids an n x n diagonal)
    b = matmul_gpu(add_vec_mat_gpu(-gmean, Z), Wsub.wt)
    alpha = matmul_gpu(b, invert_mat(crossprod_gpu(Wsub.wt, Wsub)))
    
    if (is.spanorm) {
      # set first column of alpha to be the same for all genes (see SpaNorm Model specification)
      a1_mean <- mean(as.matrix(alpha)[, 1])
      if (checkGPU() && is_torch_tensor(alpha) && backend %in% c("gpu", "auto")) {
        # Build constant first column on GPU (inherit alpha's exact dtype/device
        # so the concat below is type-consistent)
        n_genes <- as.integer(dim(alpha)[[1]])
        alpha_first <- torch::torch_multiply(
          torch::torch_ones(c(n_genes, 1L), dtype = alpha$dtype, device = alpha$device),
          torch::torch_tensor(as.numeric(a1_mean), dtype = alpha$dtype, device = alpha$device)
        )

        # Slice first and remaining columns of W on GPU (narrow preserves the 2D
        # rank; plain `[` indexing may drop the singleton column dimension)
        n_cov <- as.integer(dim(Wsub)[[2]])
        W_first <- torch::torch_narrow(Wsub, dim = 2L, start = 1L, length = 1L)

        Wa1 <- tcrossprod_gpu(alpha_first, W_first)

        W_rest <- torch::torch_narrow(Wsub, dim = 2L, start = 2L, length = n_cov - 1L)

        W_rest.wt <- mult_vec_mat_gpu(wt.cell, W_rest) # row-scale W by cell weights (avoids an n x n diagonal)
        b <- matmul_gpu(add_vec_mat_gpu(-gmean, Z) - Wa1, W_rest.wt)
        alpha_other <- matmul_gpu(b, invert_mat(crossprod_gpu(W_rest.wt, W_rest) + lambda.a))

        alpha <- torch::torch_cat(list(alpha_first, alpha_other), dim = 2L)
        rm(Wa1)
      } else {
        # CPU/base fallback
        alpha[, 1] <- rep(a1_mean, nrow(alpha))
        Wa1 <- tcrossprod_gpu(
          toGPUMatrix(matrix(alpha[, 1], ncol = 1), backend = backend),
          toGPUMatrix(matrix(Wsub[, 1], ncol = 1), backend = backend)
        )

        W_rest <- Wsub[, -1, drop = FALSE]
        W_rest.wt <- mult_vec_mat_gpu(wt.cell, W_rest) # row-scale W by cell weights (avoids an n x n diagonal)
        b <- matmul_gpu(add_vec_mat_gpu(-gmean, Z) - Wa1, W_rest.wt)
        alpha_other <- matmul_gpu(b, invert_mat(crossprod_gpu(W_rest.wt, W_rest) + lambda.a))
        alpha <- cbind(alpha[, 1, drop = FALSE], as.matrix(alpha_other))
        rm(Wa1)
      }
    }

    # if new alpha est has missing/inf values, revert to previous estimate
    # (scalar reduction on-device for tensors; no full materialisation)
    if (hasBadValues(alpha)) alpha <- alpha.old

    # reduce outliers: winsorise columns to median +/- 4*MAD. winsoriseCols keeps
    # tensors on-device (avoids the GPU->CPU->GPU round trip) while reproducing
    # the exact matrixStats behaviour on the CPU path.
    alpha = winsoriseCols(alpha, k = 4)

    # step2: update gmean
    Z.res = Z - tcrossprod_gpu(alpha, Wsub)
    gmean = rowSums_gpu(Z.res * sig.inv) / rowSums_gpu(sig.inv)

    # calculate current logl
    lmu.hat = add_vec_mat_gpu(gmean, tcrossprod_gpu(alpha, Wsub), backend = backend)
    loglik.tmp = as.numeric(colSums_gpu(dnbinom_gpu(Ysub, mu = exp(lmu.hat), size = 1 / psi, log = TRUE)))

    # check degenerate case
    degener = sum(loglik) > sum(loglik.tmp)
    if (is.na(degener) | degener) {
      # alter stepsize of genes
      check.gene = loglik > loglik.tmp
      step[check.gene] = step[check.gene] * step.factor
      # revert to previous best
      gmean = best.gmean
      alpha = best.a
      # increase halving counter
      halving = halving + 1
      if (halving >= 3) {
        loglik.tmp = loglik
        degener = !degener
      }
    }

    # check degener (again in case halving exceeds maximum, i.e., 3)
    if (!degener) {
      loglik = loglik.tmp
      logl.beta = c(sum(loglik), logl.beta) # push to stack
      iter = iter + 1
      halving = 0 # reset
    }

    # check convergence
    conv.logl = FALSE
    if (iter > 2) { # should be 2 because iter is post-incremented
      conv.logl = ifelse((logl.beta[1] - logl.beta[2]) / abs(logl.beta[2]) < tol, TRUE, FALSE)
    }
    conv = conv.logl | iter > maxit.nb
  }

  # convergence message
  if (iter > maxit.nb) {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (did not converged)", iter, logl.beta[1]))
  } else {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (converged)", iter, logl.beta[1]))
  }

  # final values
  fit.nb = list(
    gmean = as.numeric(gmean),
    alpha = toRMatrix(alpha),
    loglik = loglik,
    loglik.iter = logl.beta
  )

  return(fit.nb)
}

normaliseLogPAC <- function(Y, scale.factor, fit.spanorm) {
  # calculate mu with all effects
  mu <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # calculate mu with biology only
  isbio <- fit.spanorm$wtype %in% "biology"
  mu.2 <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])
  # winsorize dispersion parameters (large disp slow the process)
  psi <- winsorisePsi(fit.spanorm$psi)

  # logPAC
  lb <- pnbinom(as.matrix(Y) - 1, mu = mu, size = 1 / psi)
  ub <- dnbinom(as.matrix(Y), mu = mu, size = 1 / psi) + lb
  p <- (lb + ub) / 2
  p <- pmax(pmin(p, 0.999), 0.001)

  # return logPAC
  normmat <- log2(qnbinom(p, mu = scale.factor * mu.2, size = 1 / psi) + 1)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

normaliseMeanBio <- function(Y, scale.factor, fit.spanorm) {
  isbio <- fit.spanorm$wtype %in% "biology"
  # calculate mu without library size effect
  normmat <- log2(calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])) # log2(mu.2)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

normaliseMeanBatch <- function(Y, scale.factor, fit.spanorm) {
  isbio = fit.spanorm$wtype %in% "batch"
  # calculate mu without library size effect
  normmat = log2(calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])) # log2(mu.2)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

normaliseMeanLS <- function(Y, scale.factor, fit.spanorm) {
  isbio = fit.spanorm$wtype %in% "ls"

  # modify W to compute effect with the median library size
  W = fit.spanorm$W[, isbio, drop = FALSE]
  W = W / W[, 1] * median(W[, 1])

  # calculate mu without library size effect
  normmat = log2(calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], W)) # log2(mu.2)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

normaliseMedianBio <- function(Y, scale.factor, fit.spanorm) {
  # calculate mu without library size effect
  isbio <- fit.spanorm$wtype %in% "biology"
  mu.2 <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])
  # winsorize dispersion parameters (large disp slow the process)
  psi <- winsorisePsi(fit.spanorm$psi)

  # median Bio
  normmat <- log2(qnbinom_gpu(0.5, mu = scale.factor * mu.2, size = 1 / psi) + 1)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

normalisePearson <- function(Y, scale.factor, fit.spanorm) {
  isbio <- fit.spanorm$wtype %in% "biology"
  # calculate mu
  mu <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # winsorize dispersion parameters (large disp slow the process)
  psi <- winsorisePsi(fit.spanorm$psi)

  # Pearson
  normmat <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, !isbio, drop = FALSE], fit.spanorm$W[, !isbio, drop = FALSE])
  normmat <- (Y - normmat) / sqrt(mu + mu^2 * psi)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

# Winsorise per-gene dispersions to exp(median(log psi) + 3*MAD), which bounds
# runtime for genes with very large dispersion. The threshold is a global
# statistic across all genes: when normalising in gene-blocks it is computed
# once on the full fit and carried as a 'psi.max' attribute on `psi` (set by
# subsetFitGenes) so block results are identical to the whole-matrix path;
# direct callers pass a plain vector and it is computed locally as before.
winsorisePsi <- function(psi) {
  psi.max <- attr(psi, "psi.max")
  if (is.null(psi.max)) {
    psi.max <- exp(median(log(psi)) + 3 * mad(log(psi)))
  }
  pmin(as.numeric(psi), psi.max)
}

calculateMu <- function(gmean, alpha, W) {
  # calculate mu (rather log of mu)
  # mu <- add_vec_mat_gpu(gmean, tcrossprod_gpu(alpha, W)) # log(mu)
  mu <- gmean + tcrossprod(alpha, W) # log(mu)
  # winsorise
  lmu.max <- matrixStats::rowMedians(mu) + 4 * matrixStats::rowMads(mu)
  mu <- exp(pmin(mu, lmu.max))

  return(mu)
}
