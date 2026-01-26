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
  Ysub.psi = as.matrix(Ysub)[, psi.idx, drop = FALSE]
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
    offs.psi = as.matrix(Wa)[, psi.idx, drop = FALSE]
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
    loglik = colSums_gpu(dnbinom_gpu(Ysub, mu = exp(add_vec_mat_gpu(gmean, Wa, backend = backend)), size = 1/psi, log = TRUE))
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
    sig.inv = 1 / (add_vec_mat_gpu(outer(psi, rep(1, nsub), "/"), exp(-lmu.hat), backend = backend))
    loglik = colSums_gpu(dnbinom_gpu(Ysub, mu = exp(lmu.hat), size = 1 / psi, log = TRUE))
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
    Z = lmu.hat + sweep(((Ysub + 0.01) / (exp(lmu.hat) + 0.01) - 1), 2, step, "*")

    # store current alpha est
    alpha.old = copy(alpha)

    # update alpha for all genes
    sig.inv = 1 / add_vec_mat_gpu(psi, exp(-lmu.hat), backend = backend)
    wt.cell = colMeans_gpu(sig.inv)
    # prevent outlier large weight
    wt.cell = toGPUVector(pmin(as.vector(wt.cell), quantile(as.vector(wt.cell), probs = 0.98)), backend = backend)
    b = add_vec_mat_gpu(-gmean, Z) %*% diag_mat(wt.cell, backend = backend) %*% Wsub
    alpha = b %*% invert_mat(crossprod_gpu(mult_vec_mat_gpu(wt.cell, Wsub), Wsub))
    
    if (is.spanorm) {
      # set first column of alpha to be the same for all genes (see SpaNorm Model specification)
      a1_mean <- mean(as.matrix(alpha)[, 1])
      if (checkGPU() && is_tf_tensor(alpha) && backend %in% c("gpu", "auto")) {
        # Build constant first column on GPU
        a_shape <- alpha$shape$as_list()
        n_genes <- as.integer(a_shape[[1]])
        dtype <- tryCatch(alpha$dtype, error = function(e) tensorflow::tf$float32)
        alpha_first <- tensorflow::tf$math$multiply(
          tensorflow::tf$ones(list(n_genes, 1L), dtype = dtype),
          tensorflow::tf$constant(as.numeric(a1_mean), dtype = dtype)
        )

        # Slice first and remaining columns of W on GPU
        w_shape <- Wsub$shape$as_list()
        n_cells <- as.integer(w_shape[[1]])
        n_cov <- as.integer(w_shape[[2]])
        W_first <- tensorflow::tf$strided_slice(Wsub, begin = list(0L, 0L), end = list(n_cells, 1L), strides = list(1L, 1L))

        Wa1 <- tcrossprod_gpu(alpha_first, W_first)

        W_rest <- tensorflow::tf$strided_slice(Wsub, begin = list(0L, 1L), end = list(n_cells, n_cov), strides = list(1L, 1L))

        b <- (add_vec_mat_gpu(-gmean, Z) - Wa1) %*% diag_mat(wt.cell, backend = backend) %*% W_rest
        alpha_other <- b %*% invert_mat(crossprod_gpu(mult_vec_mat_gpu(wt.cell, W_rest), W_rest) + lambda.a)

        alpha <- tensorflow::tf$concat(list(alpha_first, alpha_other), axis = 1L)
        rm(Wa1)
      } else {
        # CPU/base fallback
        alpha[, 1] <- rep(a1_mean, nrow(alpha))
        Wa1 <- tcrossprod_gpu(
          toGPUMatrix(matrix(alpha[, 1], ncol = 1), backend = backend),
          toGPUMatrix(matrix(Wsub[, 1], ncol = 1), backend = backend)
        )

        W_rest <- Wsub[, -1, drop = FALSE]
        b <- (add_vec_mat_gpu(-gmean, Z) - Wa1) %*% diag_mat(wt.cell, backend = backend) %*% W_rest
        alpha_other <- b %*% invert_mat(crossprod_gpu(mult_vec_mat_gpu(wt.cell, W_rest), W_rest) + lambda.a)
        alpha <- cbind(alpha[, 1, drop = FALSE], as.matrix(alpha_other))
        rm(Wa1)
      }
    }

    # if new alpha est has missing/inf values, revert to previous estimate
    tmp = as.matrix(alpha)
    if (any(is.na(tmp)) || any(is.nan(tmp)) || any(is.infinite(tmp))) alpha <- alpha.old
    rm(tmp)

    # reduce outliers
    alpha = as.matrix(alpha)
    a.med = matrixStats::colMedians(alpha)
    a.mad = matrixStats::colMads(alpha)
    a.max = a.med + 4 * a.mad
    a.min = a.med - 4 * a.mad
    alpha = pmin(t(alpha), a.max) # +ve outliers, orig: t(pmin(t(alpha), a.max))
    alpha = t(pmax(alpha, a.min)) # -ve outliers, orig: t(pmax(t(alpha), a.min))
    alpha = toGPUMatrix(alpha, backend = backend)

    # step2: update gmean
    Z.res = Z - tcrossprod_gpu(alpha, Wsub)
    gmean = rowSums_gpu(Z.res * sig.inv) / rowSums_gpu(sig.inv)

    # calculate current logl
    lmu.hat = add_vec_mat_gpu(gmean, tcrossprod_gpu(alpha, Wsub), backend = backend)
    loglik.tmp = colSums_gpu(dnbinom_gpu(Ysub, mu = exp(lmu.hat), size = 1 / psi, log = TRUE))

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
    gmean = as.vector(gmean),
    alpha = toRMatrix(alpha),
    loglik = loglik,
    loglik.iter = logl.beta
  )

  return(fit.nb)
}

normaliseLogPAC <- function(Y, scale.factor, fit.spanorm, chunk_size = 5000) {
  # Winsorize dispersion parameters (small, done once)
  psi <- fit.spanorm$psi
  psi.max <- exp(median(log(psi)) + 3 * mad(log(psi)))
  psi <- pmin(psi, psi.max)
  
  # Identify biology columns
  isbio <- fit.spanorm$wtype %in% "biology"
  
  # Pre-compute winsorization thresholds on full dataset for numerical consistency
  n_spots <- ncol(Y)
  
  # Compute log(mu) for full dataset to get winsorization limits
  lmu_full <- fit.spanorm$gmean + tcrossprod(fit.spanorm$alpha, fit.spanorm$W)
  lmu.max_full <- matrixStats::rowMedians(lmu_full) + 4 * matrixStats::rowMads(lmu_full)
  
  lmu_full_bio <- fit.spanorm$gmean + tcrossprod(fit.spanorm$alpha[, isbio, drop = FALSE], 
                                                  fit.spanorm$W[, isbio, drop = FALSE])
  lmu.max_bio <- matrixStats::rowMedians(lmu_full_bio) + 4 * matrixStats::rowMads(lmu_full_bio)
  
  rm(lmu_full, lmu_full_bio)
  gc(verbose = FALSE)
  
  # Process in chunks to avoid memory spike
  n_chunks <- ceiling(n_spots / chunk_size)
  
  if (n_chunks > 1) {
    message(sprintf("  Processing %d spots in %d chunks (chunk_size=%d) to manage memory", 
                    n_spots, n_chunks, chunk_size))
  }
  
  # Pre-allocate final matrix (more memory-efficient than list)
  n_genes <- nrow(Y)
  normmat <- matrix(0, nrow = n_genes, ncol = n_spots)
  
  for (i in seq_len(n_chunks)) {
    # Define chunk indices
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_spots)
    chunk_idx <- start_idx:end_idx
    
    if (n_chunks > 1 && i %% max(1, floor(n_chunks / 10)) == 0) {
      message(sprintf("    Processing chunk %d/%d (spots %d-%d)", 
                      i, n_chunks, start_idx, end_idx))
    }
    
    # Subset data for this chunk
    Y_chunk <- Y[, chunk_idx, drop = FALSE]
    W_chunk <- fit.spanorm$W[chunk_idx, , drop = FALSE]
    
    # Calculate mu with all effects (using pre-computed winsorization limit)
    lmu_chunk <- fit.spanorm$gmean + tcrossprod(fit.spanorm$alpha, W_chunk)
    lmu_chunk <- exp(pmin(lmu_chunk, lmu.max_full))  # reuse variable to save memory
    
    # Calculate mu with biology only (using pre-computed winsorization limit)
    lmu.2_chunk <- fit.spanorm$gmean + tcrossprod(fit.spanorm$alpha[, isbio, drop = FALSE], 
                                                   W_chunk[, isbio, drop = FALSE])
    lmu.2_chunk <- exp(pmin(lmu.2_chunk, lmu.max_bio))
    
    # logPAC computation (chunk only)
    Y_chunk <- as.matrix(Y_chunk)
    lb <- pnbinom(Y_chunk - 1, mu = lmu_chunk, size = 1 / psi)
    ub <- dnbinom(Y_chunk, mu = lmu_chunk, size = 1 / psi) + lb
    p <- (lb + ub) / 2
    p <- pmax(pmin(p, 0.999), 0.001)
    
    # Compute normalized values and store directly in final matrix
    normmat[, start_idx:end_idx] <- log2(qnbinom(p, mu = scale.factor * lmu.2_chunk, size = 1 / psi) + 1)
    
    # Clean up chunk-specific objects
    rm(lmu_chunk, lmu.2_chunk, lb, ub, p, Y_chunk, W_chunk)
    gc(verbose = FALSE)
  }
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
  psi <- fit.spanorm$psi
  psi.max <- exp(median(log(psi)) + 3 * mad(log(psi)))
  psi <- pmin(psi, psi.max)

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
  psi <- fit.spanorm$psi
  psi.max <- exp(median(log(psi)) + 3 * mad(log(psi)))
  psi <- pmin(psi, psi.max)

  # Pearson
  normmat <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, !isbio, drop = FALSE], fit.spanorm$W[, !isbio, drop = FALSE])
  normmat <- (Y - normmat) / sqrt(mu + mu^2 * psi)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
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
