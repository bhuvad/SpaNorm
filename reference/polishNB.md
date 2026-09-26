# Converge every gene of a shared negative binomial fit to its own optimum

[`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md) fits
every gene in one IRLS loop: it shares a single gene-averaged cell
weight vector across genes, decides step-halving and convergence on the
aggregate log-likelihood, and clamps coefficients across genes. For most
genes the result is indistinguishable from the gene's own optimum, but a
bright, cell-type-restricted gene can be left well short of it, with an
inflated dispersion. `polishNB()` takes such a fit and converges each
gene separately, by damped Newton on that gene's own penalised negative
binomial log-likelihood \$\$\sum_i \log f\_{NB}(y_i; \mu_i, \psi) -
\frac{1}{2} \sum_j \lambda_j \alpha_j^2, \qquad \log \mu = W \alpha +
\mathrm{offset}.\$\$

## Usage

``` r
polishNB(
  Y,
  W,
  alpha,
  psi,
  lambda.a = 0,
  offset = NULL,
  absorb = NULL,
  absorb.batch = NULL,
  start.cols = NULL,
  psi.method = c("profile", "fixed"),
  psi.range = c(0.001, 1000),
  warm = FALSE,
  maxit = 50L,
  tol = 1e-08,
  engine = c("batch", "gene"),
  batch.size = NULL,
  block.size = NULL,
  backend = c("cpu", "auto", "gpu"),
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## Arguments

- Y:

  a genes x cells matrix of integer counts (dense, sparse or
  DelayedArray; densified one gene block at a time). The negative
  binomial likelihood is undefined on non-integer values, so a
  non-integer assay is refused.

- W:

  a cells x p numeric design matrix.

- alpha:

  a genes x p matrix of starting coefficients, typically
  `fitNB()$alpha`.

- psi:

  the starting per-gene dispersions (length `nrow(Y)`, or one value for
  every gene), typically `fitNB()$psi`.

- lambda.a:

  the ridge penalty, a single value or one per column of `W`. It is
  applied as given: each gene's objective subtracts
  `0.5 * sum(lambda.a * alpha^2)`, with no scaling by the number of
  cells or genes. The caller owns the scaling, so a fit made with a
  scaled penalty must pass the scaled values here.

- offset:

  `NULL` (the default), a log-scale offset with one value per cell (the
  same for every gene), or a genes x cells matrix of them, added to
  every linear predictor with its coefficient fixed at 1:
  `log(mu) = W alpha + offset`. It must be finite. A `Matrix` offset is
  made a base matrix; a torch tensor is refused (pass it on the host; on
  a device the engine moves it). A matrix is held dense, so prefer a
  vector when the offset is the same for every gene.

- absorb:

  which columns of `W` the per-gene Newton solver absorbs by a Schur
  complement (exact; it makes a wide indicator block cost one
  dense-column gram per iteration). `NULL` (the default) absorbs
  nothing. A logical over the columns of `W` makes each marked column
  its own 1x1 block; the marked columns must be 0/1 indicators that
  partition the cells. A per-column grouping (integer, factor or
  character ids, `NA` for a dense column) makes the columns sharing an
  id one block; no cell may load on two blocks. See
  [`nbNewtonSolver()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md).

- absorb.batch:

  the absorption for the shared-factor batched solver, which is what
  runs on a device (`backend = "gpu"`, or `"auto"` with a GPU found):
  `NULL` or a logical over the columns of `W`. That solver absorbs 1x1
  blocks only, so it cannot take a grouping with multi-column blocks.
  `NULL` (the default) passes a logical `absorb` through unchanged and
  sends a grouping to the dense batched solver (exact, only slower); a
  caller that knows the 1x1 subset of its grouping (e.g. the nested
  indicators inside a random-slope fit's per-sample blocks) passes it
  here. Not used on the CPU.

- start.cols:

  a logical over the columns of `W` marking the indicator columns (such
  as cell-type intercepts) that the sane start fills with the gene's log
  mean over that column's cells, or `NULL` to put the overall log mean
  on the first column.

- psi.method:

  how the dispersion is set at the converged mean: `"profile"` (profile
  maximum likelihood per gene) or `"fixed"` (the input `psi` is kept and
  only the mean is converged).

- psi.range:

  the search interval for the profile dispersion, two increasing
  positive numbers. A gene whose optimum falls on either end keeps its
  input dispersion.

- warm:

  logical; `alpha` and `psi` are an already converged fit at a nearby
  penalty. A warm polish is a few damped Newton steps at the held
  dispersion, with no dispersion search and no restart check.

- maxit, tol:

  the Newton iteration cap and the relative log-likelihood tolerance,
  per gene.

- engine:

  `"batch"` (the default) runs a batch of genes through each Newton step
  so they share every read of the design; `"gene"` is the per-gene
  reference implementation. Both reach the same optimum; the batched
  profile dispersion is a bisection and the per-gene one
  [`optimize()`](https://rdrr.io/r/stats/optimize.html), so they agree
  on `psi` to [`optimize()`](https://rdrr.io/r/stats/optimize.html)'s
  tolerance.

- batch.size:

  genes per batched Newton (`engine = "batch"`), or `NULL` to size it
  from a per-worker memory budget,
  `options(SpaNorm.polish.mem.budget = <bytes>)` (default 1e9; the older
  `spiDE.polish.mem.budget` is read as a fallback).

- block.size:

  genes per block, the unit of densification and dispatch, or `NULL` for
  at least one block per worker and at most 2,000 genes.

- backend:

  `"cpu"` (the default), `"auto"` or `"gpu"`. A device needs
  `engine = "batch"` and float64, and forces serial dispatch.

- BPPARAM:

  a `BiocParallelParam` over gene blocks. With more than one worker and
  the RhpcBLASctl package installed, each worker runs its BLAS and
  OpenMP single-threaded, because forked workers inherit the parent's
  thread count and oversubscribing the cores costs an order of magnitude
  per gene. Without RhpcBLASctl the workers keep the thread count they
  inherit; install it, or set the BLAS threads to 1 (e.g.
  `OPENBLAS_NUM_THREADS=1`) before starting R, when using several
  workers.

- verbose:

  logical; report progress.

## Value

a list with `alpha` (genes x p), `psi` and `loglik` (the penalised
log-likelihood at the returned fit, `NA` for a gene that was not
polished), and `polish`, a data frame with one row per gene (row names
from `alpha`): `iterations` (Newton steps), `psi_fitnb` (the input
dispersion), `restarted` (the sane start was used), `capped` (a Newton
pass hit `maxit`), `singular` (a singular information matrix),
`psi_bound` (the dispersion optimum was on its search bound, so the
input value was kept) and `polished` (`FALSE` when the gene kept its
input fit).

## Details

With `psi.method = "profile"` the dispersion is then re-estimated by
profile maximum likelihood at the converged mean and the mean
re-polished, twice; a gene whose dispersion optimum sits on the search
bound keeps its input dispersion (`psi_bound`). A gene whose starting
point is degenerate (a fitted log-mean below -10 at a cell with a
positive count) or whose Newton diverges restarts from a sane point: the
log mean over the cells of each `start.cols` column, or the overall log
mean on the first column. A gene that cannot be polished from either
start keeps its input coefficients and dispersion, with
`polished = FALSE`.

Once the shared fit is made the genes are independent, so this stage is
exact when blocked: the counts are densified one gene block at a time
(`block.size`) and the blocks are dispatched with `BPPARAM`. The
dispersion moderation in
[`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md) works
across genes, which is why that fit must see the whole gene set and this
one need not.

The integer-count check reads only the first 20 genes. It is a cheap
guard against an assay that is non-integer throughout (a back-transform
such as `2^logcounts - 1`, on which every gene's dispersion would
silently run to its upper bound), not a scan of every value: a
non-integer count in a later gene is not detected.

## See also

[`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md),
[`nbNewtonSolver()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md).

## Examples

``` r
set.seed(1)
W <- cbind(1, rnorm(200))
Y <- t(replicate(5, rnbinom(200, mu = exp(1 + 0.3 * W[, 2]), size = 3)))
fit <- fitNB(Y, W, verbose = FALSE, backend = "cpu")
pol <- polishNB(Y, W, fit$alpha, fit$psi)
pol$polish
#>   iterations psi_fitnb restarted capped singular psi_bound polished
#> 1          5 0.3049685     FALSE  FALSE    FALSE     FALSE     TRUE
#> 2          5 0.3522451     FALSE  FALSE    FALSE     FALSE     TRUE
#> 3          5 0.3532502     FALSE  FALSE    FALSE     FALSE     TRUE
#> 4          4 0.2985352     FALSE  FALSE    FALSE     FALSE     TRUE
#> 5          4 0.3022655     FALSE  FALSE    FALSE     FALSE     TRUE
```
