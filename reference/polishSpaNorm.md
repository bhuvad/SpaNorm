# Converge each gene of a SpaNorm fit to its own optimum

[`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
fits every gene in one shared IRLS loop, with a gene-averaged cell
weight vector and an aggregate convergence criterion, so a gene can be
left short of its own optimum (bright, spatially restricted genes most
of all). `polishSpaNorm()` takes the stored fit and converges each gene
separately under SpaNorm's own model, stores the polished fit, and
rewrites the normalised assay from it.

## Usage

``` r
polishSpaNorm(
  spe,
  adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
  scale.factor = 1,
  psi.method = c("fixed", "profile"),
  ls = c("fixed", "joint"),
  cells = c("all", "fit"),
  null = TRUE,
  maxit = 50L,
  tol = 1e-08,
  engine = c("batch", "gene"),
  batch.size = NULL,
  backend = c("cpu", "auto", "gpu"),
  BPPARAM = BiocParallel::SerialParam(),
  overwrite = FALSE,
  verbose = TRUE,
  assay = NULL,
  ...
)

# S4 method for class 'SpatialExperiment'
polishSpaNorm(
  spe,
  adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
  scale.factor = 1,
  psi.method = c("fixed", "profile"),
  ls = c("fixed", "joint"),
  cells = c("all", "fit"),
  null = TRUE,
  maxit = 50L,
  tol = 1e-08,
  engine = c("batch", "gene"),
  batch.size = NULL,
  backend = c("cpu", "auto", "gpu"),
  BPPARAM = BiocParallel::SerialParam(),
  overwrite = FALSE,
  verbose = TRUE,
  assay = NULL,
  ...
)

# S4 method for class 'Seurat'
polishSpaNorm(
  spe,
  adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
  scale.factor = 1,
  psi.method = c("fixed", "profile"),
  ls = c("fixed", "joint"),
  cells = c("all", "fit"),
  null = TRUE,
  maxit = 50L,
  tol = 1e-08,
  engine = c("batch", "gene"),
  batch.size = NULL,
  backend = c("cpu", "auto", "gpu"),
  BPPARAM = BiocParallel::SerialParam(),
  overwrite = FALSE,
  verbose = TRUE,
  assay = NULL,
  ...
)
```

## Arguments

- spe:

  a SpatialExperiment or Seurat object that
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
  has been run on, with the raw counts in its 'counts' assay (layer, for
  Seurat).

- adj.method:

  a character, specifying the method used to rewrite the normalised data
  from the polished fit (default 'auto'), as in
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md).

- scale.factor:

  a numeric, specifying the scaling factor for the adjusted counts, as
  in
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md).

- psi.method:

  how each gene's dispersion is set: `"fixed"` (the default) keeps the
  fit's dispersion and converges the mean under it; `"profile"`
  re-estimates it by profile maximum likelihood at the converged mean
  and re-polishes (see
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)).

- ls:

  the library-size coefficient shared by every gene: `"fixed"` (the
  default) holds it at the fit's value; `"joint"` re-estimates it
  jointly with the per-gene coefficients, at the optimum of the total
  penalised likelihood over the polished genes (see Details).

- cells:

  the cells/spots each gene is converged over: `"all"` (the default), or
  `"fit"` for the ones
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
  sampled to fit the model.

- null:

  a logical, specifying whether to also polish the null model stored by
  [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  ('SpaNormNull'), when there is one (default TRUE). The SVG test
  compares the two fits, so they must be estimated the same way.

- maxit, tol:

  the Newton iteration cap and the relative log-likelihood tolerance,
  per gene.

- engine, batch.size, backend:

  passed to
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md):
  the batched or per-gene Newton engine, genes per batch, and the
  compute backend.

- BPPARAM:

  a BiocParallelParam object over gene blocks, for the polish (see
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md))
  and the normalisation (default
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)).

- overwrite:

  a logical. A fit that is already polished is refused unless this is
  TRUE (default FALSE); with TRUE it is polished again, starting from
  the unpolished fit kept as 'SpaNormUnpolished' (and
  'SpaNormNullUnpolished'), never from the polished one.

- verbose:

  a logical, specifying whether to show update messages (default TRUE).

- assay:

  a character, specifying the assay holding the raw counts for Seurat
  objects (default NULL uses the object's default assay). Ignored for
  SpatialExperiment objects.

- ...:

  further arguments passed to
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md),
  such as `block.size` or `psi.range`.

## Value

a SpatialExperiment or Seurat object holding the polished fit(s), with
the normalised data rewritten in 'logcounts' or 'data', respectively.

## Details

SpaNorm's mean for gene \\g\\ in cell/spot \\i\\ is \$\$\log \mu\_{gi} =
\bar{\mu}\_g + a_1 w\_{i1} + \sum\_{j \ge 2} W\_{ij} \alpha\_{gj},\$\$
where \\w\_{i1}\\ is the log library size and its coefficient \\a_1\\ is
shared by every gene. The polish maximises each gene's penalised
negative binomial log-likelihood by damped Newton (see
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)),
with \\a_1\\ held at the fit's value (`ls = "fixed"`), \\\bar{\mu}\_g\\
a per-gene unpenalised intercept, and the other coefficients
ridge-penalised as
[`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
penalises them: by `lambda.a` for the biology and library-size splines,
not at all for batch, times the number of cells/spots.

The polished coefficients are the optimum of the unwinsorised penalised
likelihood. The normalisation (and
[`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md))
still applies SpaNorm's usual winsorisation, exactly as for an
unpolished fit: it forms the mean through
[`calculateMu()`](https://bhuvad.github.io/spaNorm/reference/calculateMu.md),
which by default winsorises the upper tail of each gene's log mean, and
caps the dispersions at `exp(median(log psi) + 4 MAD(log psi))`.

A gene with no counts in the polished cells/spots has no finite optimum
(its intercept runs to minus infinity), and nor has a gene with no
counts in every polished cell/spot of one batch level (that level's
unpenalised coefficient runs to minus infinity), so neither is polished:
it keeps its input coefficients and dispersion, and its diagnostics row
reads `polished = FALSE` with `iterations = 0` and `held_out` naming the
reason, `"all-zero"` or `"zero-batch-level"` (`NA` for a gene that was
polished). Genes held out for a batch level are counted in a warning.
The levels are the groups of polished cells/spots that share one row of
the model's unpenalised columns, the intercept and the batch indicators:
for a batch factor, its levels; for several batch variables, each
combination of levels present, which can hold out a gene whose optimum
under an additive batch design is finite. Groups are formed only when
every unpenalised column other than the intercept is a 0/1 indicator;
when one is not (a numeric batch covariate, or a spline left unpenalised
by `lambda.a = 0`), only the all-zero rule applies. A gene the Newton
engine cannot converge from either start also keeps its input fit (see
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md));
a warning counts these genes, and how many had a singular information
matrix (for example, a batch matrix given with every level's indicator,
so that the batch columns are collinear with the intercept). The fit is
marked polished either way.

With `ls = "joint"` the shared coefficient \\a_1\\ is moved to the joint
optimum of the total penalised log-likelihood over the polished genes:
after the per-gene polish at the fit's \\a_1\\, Newton steps on the
profile log-likelihood of \\a_1\\ (each gene's own coefficients profiled
out; the information is the Fisher information of the polish) alternate
with a warm re-polish of every gene at the candidate value, and a step
is halved until the total does not fall. It stops when the standardised
score \\\|U\|/\sqrt{I}\\ (a1's distance from its optimum in units of its
own profiled SE) is below `1e-3`, or after 10 steps; `maxit` and `tol`
are the per-gene Newton's, in the cold pass and every re-polish alike.
With `psi.method = "profile"` each re-polish also re-profiles the
dispersion. A gene whose information is singular at a step is left out
of that step's score, information and total. If the profiled information
is not positive (\\a_1\\ not identified: the log library size is
collinear with the model's other terms), the fit is the `ls = "fixed"`
polish, with a warning. The joint estimate weights genes by their
information, so bright genes dominate it, where
[`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)'s
\\a_1\\ is an unweighted mean over genes. Genes that are not polished
(held out, or not convergeable) do not inform \\a_1\\; they keep their
input \\\bar{\mu}\_g\\, other coefficients and dispersion exactly, but
take the new shared \\a_1\\, since it is one value for every gene, so
their fitted mean shifts by \\(a_1^{new} - a_1^{old}) w\_{i1}\\. With no
polished gene at all, \\a_1\\ is left at the fit's value with a warning.

The polished fit replaces 'SpaNorm' in the object's metadata (`@misc`
for Seurat) and the input fit is kept as 'SpaNormUnpolished'. With
`null = TRUE` a stored 'SpaNormNull' is polished with the same settings
and the input kept as 'SpaNormNullUnpolished'. SVG results in `rowData`
were computed from the unpolished fit, so they are removed with a
warning; rerun
[`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md).
The normalised assay ('logcounts', or the 'data' layer for Seurat) is
rewritten from the polished fit.

The polished fit's `polish` slot (see
[`isPolished()`](https://bhuvad.github.io/spaNorm/reference/isPolished.md))
holds `settings` (`psi.method`, `ls`, `cells`, `maxit`, `tol`, the
penalty vector `pen`, the library-size coefficient before (`a1.input`)
and after (`a1`), with `psi.method = "profile"` the dispersion's search
interval `psi.range` (as passed through `...`, else
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)'s
default), and the SpaNorm version) and `genes`, one row per gene:
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)'s
diagnostics plus `loglik`, the penalised log-likelihood over the
polished cells at the returned fit (`NA` for a gene that was not
polished), and `held_out` (see above). The fit's `loglik` slot is left
as the shared fit's iteration trace. With `ls = "joint"`, `settings`
also holds `ls.iterations` (accepted steps on \\a_1\\), `ls.maxit` and
`ls.tol` (the cap on those steps and the standardised-score stop, 10 and
`1e-3`, the latter in units of a1's profiled SE), `ls.score` (the final
pooled score \\U\\), `ls.se` (\\1/\sqrt{I}\\, the profiled standard
error of \\a_1\\), `ls.singular` (genes left out of \\U\\ and \\I\\ at
some step because their information was singular) and `ls.converged`,
and the per-gene `iterations`, `capped` and `singular` include the warm
re-polishes.

The negative binomial likelihood is defined on counts only, so a
non-integer counts assay (for example a back-transform such as
`2^logcounts - 1`) is refused.

## See also

[`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md),
[`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md),
[`isPolished()`](https://bhuvad.github.io/spaNorm/reference/isPolished.md).

## Examples

``` r
data(HumanDLPFC)
# \donttest{
top <- order(-Matrix::rowSums(SummarizedExperiment::assay(HumanDLPFC, "counts")))[1:50]
spe <- SpaNorm(HumanDLPFC[top, ], sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> (1/2) Fitting SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -28161.783816
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -28161.783816
#> iter:  1, iter:  2, log-likelihood: -27103.362846
#> iter:  1, iter:  3, log-likelihood: -27022.731294
#> iter:  1, iter:  4, log-likelihood: -27015.689392
#> iter:  1, iter:  4, log-likelihood: -27015.689392
#> iter:  1, iter:  4, log-likelihood: -27015.689392
#> iter:  1, iter:  5, log-likelihood: -27015.689392 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -26999.708156
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -26999.708156
#> iter:  2, iter:  1, log-likelihood: -26999.708156
#> iter:  2, iter:  1, log-likelihood: -26999.708156
#> iter:  2, iter:  2, log-likelihood: -26999.708156
#> iter:  2, iter:  2, log-likelihood: -26999.708156
#> iter:  2, iter:  2, log-likelihood: -26999.708156
#> iter:  2, iter:  3, log-likelihood: -26999.708156 (converged)
#> iter:  3, log-likelihood: -26999.708156 (converged)
#> (2/2) Normalising data
spe <- polishSpaNorm(spe)
#> (1/2) Polishing SpaNorm model
#>   converging 50 genes in batches of 4447 (1 block)
#>     block 1/1 (0.0 min elapsed)
#>   polished 50/50 genes (0 restarted, 0 not converged, 0 singular, 0 dispersion at a bound)
#> (2/2) Normalising data
isPolished(S4Vectors::metadata(spe)$SpaNorm)
#> [1] TRUE
# }
```
