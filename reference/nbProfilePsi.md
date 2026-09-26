# Profile-ML dispersion at fixed coefficients, blocked over genes

For each gene, the negative binomial dispersion that maximises the
likelihood at the mean `exp(W alpha + offset)`, the coefficients held
fixed: one profile search per gene and no Newton, by the same bisection
on the score that
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)'s
batched engine uses. A gene whose optimum sits on a bound of `psi.range`
(an under-dispersed or near-empty gene), or whose search is not finite,
keeps its input dispersion rather than a boundary value stored as an
estimate, the rule
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)
applies. A caller that re-polishes the mean at a held dispersion
(`polishNB(warm = TRUE)`) uses this to put the dispersion back at the
reported mean.

## Usage

``` r
nbProfilePsi(
  Y,
  W,
  alpha,
  psi,
  psi.range = c(0.001, 1000),
  block.size = NULL,
  BPPARAM = BiocParallel::SerialParam(),
  offset = NULL
)
```

## Arguments

- Y:

  a genes x cells matrix of counts (dense, sparse or DelayedArray;
  densified one gene block at a time).

- W:

  a cells x p numeric design matrix.

- alpha:

  a genes x p matrix of coefficients.

- psi:

  the per-gene dispersions (length `nrow(Y)`, or one value for every
  gene), kept for a gene whose optimum is on a bound.

- psi.range:

  the search interval for the dispersion.

- block.size, BPPARAM:

  gene blocking and dispatch, as in
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md).

- offset:

  `NULL`, a per-cell log-scale offset (length `ncol(Y)`) or a genes x
  cells matrix, added to every linear predictor. The dispersion is
  profiled at the mean it gives, so the offset must be the one the
  coefficients were fitted with.

## Value

a numeric vector of dispersions, one per gene.

## See also

[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md).

## Examples

``` r
set.seed(1)
W <- cbind(1, rnorm(200))
Y <- t(replicate(3, rnbinom(200, mu = exp(1 + 0.3 * W[, 2]), size = 4)))
nbProfilePsi(Y, W, matrix(c(1, 0.3), 3, 2, byrow = TRUE), rep(0.5, 3))
#> [1] 0.2917504 0.2388419 0.2547512
```
