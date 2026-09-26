# Penalised NB Newton solvers and batched grams (for package developers)

Low-level building blocks of
[`polishNB`](https://bhuvad.github.io/spaNorm/reference/polishNB.md),
exported so that a downstream package (e.g. spiDE) can form the
penalised covariance of a polished fit from the SAME factorisation the
polish used. The information is `crossprod(W, w * W) + diag(pen)`.
`absorb` marks an indicator block whose part of the information is
block-diagonal and is eliminated by a Schur complement: `NULL` (dense),
a logical over `W`'s columns (each marked column its own 1x1 block), or
an integer block id per column (`NA` = dense).

## Usage

``` r
nbNewtonSolver(W, pen, absorb = NULL)

nbNewtonSolverBatch(W, pen, absorb = NULL)

nbGramBatch(
  W,
  wt_block,
  penalty_diag = NULL,
  backend = "cpu",
  cell.tile = NULL
)

nbAbsorbGramBatch(W, pen, absorb, wt_block, cell.tile = NULL, parts = FALSE)
```

## Arguments

- W:

  a cells x p design (base matrix, or torch tensor for the batch forms).

- pen:

  a length-p ridge penalty.

- absorb:

  see Description.

- wt_block:

  a genes x cells weight matrix (one row per gene).

- penalty_diag:

  `NULL`, or a length-p ridge penalty added to the diagonal of every
  gene's gram.

- backend:

  the resolved backend; unused on base R matrices, kept for symmetry
  with the other batched helpers.

- cell.tile:

  cells per accumulation tile, or `NULL` for all at once. A gram is a
  sum over cells, so every tiling gives the same result; a tile bounds
  the `batch x cells x p` weighted design on a device.

- parts:

  if `TRUE`, return the pieces a Newton step's back-substitution needs
  (`S`, `B`, `cvec`, `xi`, `zi`) rather than the Schur complements `S`
  alone.

## Value

`nbNewtonSolver()`: a list of closures `factor(w)`, `solve(state, s)`
(the Newton step, `NULL` if singular) and `xcov(state)` (the dense-block
covariance). The batch forms return the batched grams or Schur
complements.

## Examples

``` r
W <- cbind(1, seq(-1, 1, length.out = 20))
s <- nbNewtonSolver(W, pen = c(0, 1))
st <- s$factor(rep(1, 20))
s$solve(st, c(1, 0))
#> [1]  5.000000e-02 -1.326682e-18
```
