# A design with a per-group block, in compact form

Describes the design matrix \$\$W = \[\\ X \mid Z_1 \mid Z_2 \mid \dots
\mid Z_G \\\],\$\$ where `X` holds the ordinary (dense) columns and the
per-group blocks \\Z_g\\ are the rows of `Z` belonging to group \\g\\
and zero elsewhere: every group has its own copy of the `q` columns of
`Z`. This is, for example, a per-patient intercept and per-patient
spatial library-size spline,
`Z = cbind(1, l * tpsBasis(x, y, df = c(3, 3)))`, with `block` the
patient. The linear predictor of a gene with coefficients `alpha` is
`X %*% alpha_x + rowSums(Z * A[block, ])`, with `A` the `G x q` matrix
of that gene's block coefficients.

## Usage

``` r
nbBlockDesign(X, Z, block, rank.tol = 1e-10)

# S3 method for class 'nbBlockDesign'
dim(x)

# S3 method for class 'nbBlockDesign'
dimnames(x)

# S3 method for class 'nbBlockDesign'
as.matrix(x, ...)

# S3 method for class 'nbBlockDesign'
print(x, ...)
```

## Arguments

- X:

  a numeric cells x `p_x` matrix of the non-absorbed columns
  (`p_x >= 1`).

- Z:

  a numeric cells x `q` matrix: each cell's values of its own group's
  block columns.

- block:

  the group of each cell, a factor (or a vector coerced with
  [`factor()`](https://rdrr.io/r/base/factor.html)) with no missing
  values. Every level is a group, including a level with no cells, whose
  coefficients are then identified by the penalty alone (or held at
  their start, see Degenerate blocks).

- rank.tol:

  the tolerance below which a group's block is treated as rank-deficient
  (see Degenerate blocks).

- x:

  an `nbBlockDesign`.

- ...:

  unused.

## Value

an object of class `nbBlockDesign`: a list holding `X`, `Z`, `block` and
precomputed indexing, with [`dim()`](https://rdrr.io/r/base/dim.html)
and [`dimnames()`](https://rdrr.io/r/base/dimnames.html) methods that
describe the implied design (so
[`nrow()`](https://rdrr.io/r/base/nrow.html),
[`ncol()`](https://rdrr.io/r/base/nrow.html) and
[`colnames()`](https://rdrr.io/r/base/colnames.html) work) and an
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) method that builds
it densely (for checking, or for a small problem).

## Details

The object is accepted as the design `W` by
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md),
[`nbProfilePsi()`](https://bhuvad.github.io/spaNorm/reference/nbProfilePsi.md)
and
[`nbNewtonSolver()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md),
which then never build the dense `n x (p_x + G q)` matrix. The `G`
blocks are absorbed exactly by a Schur complement, as
`polishNB(absorb = )` absorbs a per-column grouping of a dense design,
so a Newton step costs one `n x (p_x + q)` gram whatever `G` is. The
results are those of the dense design `as.matrix(W)` with the block
columns absorbed by group
(`absorb = c(rep(NA, p_x), rep(seq_len(G), each = q))`).

**Coefficient layout.** The implied design has `p = p_x + G * q`
columns, in the order of
[`as.matrix()`](https://rdrr.io/r/base/matrix.html): the columns of `X`,
then group 1's `q` columns, then group 2's, and so on (groups in the
order of `levels(block)`). `alpha`, `lambda.a` and `start.cols` follow
that layout; the per-group coefficients of gene `g` are
`matrix(alpha[g, -(1:p_x)], G, q, byrow = TRUE)`. `lambda.a` and
`start.cols` may also be given for `[X | Z]` only (length `p_x + q`), in
which case the `Z` part is used for every group.

**Degenerate blocks.** Each group's block of the information matrix,
\\C_g = Z_g' \mathrm{diag}(w) Z_g + \mathrm{diag}(\lambda_g)\\, is
scaled to unit diagonal and Cholesky-factorised. A group whose block is
rank-deficient or nearly so – fewer cells than columns, a column
constant or zero over the group's cells, a level of `block` with no
cells at all, any of these with a zero penalty – is recognised by its
smallest squared Cholesky pivot falling below `rank.tol` (on the
unit-diagonal scale that pivot is \\1 - R^2\\ of a column regressed on
the columns before it, so the default `1e-10` flags a variance inflation
above \\10^{10}\\). Such a block is inverted by a truncated
eigendecomposition instead (a generalised inverse; eigenvalues below
`rank.tol` times the largest are dropped), and a column with no weight
and no penalty in the group is dropped outright. The directions this
drops are combinations of the group's coefficients that do not change
the linear predictor on its cells and carry no penalty, so the objective
is flat along them: the Newton step is the minimum-norm step in the
unit-diagonal scaling, which leaves the coefficients along those
directions at their starting values (an unpenalised column of a group
with no cells keeps its start exactly). The likelihood, the fitted
means, the dispersion, the dense-column coefficients and their
covariance are those of the same design with the aliased columns
removed; the individual block coefficients of an aliased group are not
identified and are one solution among many. The dense grouped path has
no such rule and reports a gene with an unpenalised aliased block as
singular.

## See also

[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md),
[`nbNewtonSolver()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md),
[`tpsBasis()`](https://bhuvad.github.io/spaNorm/reference/tpsBasis.md).

## Examples

``` r
set.seed(1)
n <- 300
patient <- factor(sample(c("P1", "P2", "P3"), n, replace = TRUE))
X <- cbind(niche = rnorm(n))
l <- rnorm(n)                                  # centred log library size
Z <- cbind(1, l)                               # per-patient intercept and slope
D <- nbBlockDesign(X, Z, patient)
dim(D)                                         # 300 x (1 + 3 * 2)
#> [1] 300   7
mu <- exp(1 + 0.3 * X[, 1] + 0.5 * l)
Y <- t(replicate(4, rnbinom(n, mu = mu, size = 5)))
a0 <- matrix(0, 4, ncol(D))
fit <- polishNB(Y, D, a0, rep(0.2, 4), lambda.a = c(0, 1e-3, 1e-3),
                start.cols = c(FALSE, TRUE, FALSE))
fit$alpha[, 1]                                 # the niche coefficient
#> [1] 0.2506974 0.2657237 0.3184350 0.3017064
```
