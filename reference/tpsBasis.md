# Thin-plate-style tensor spline basis over spatial coordinates

The spatial basis SpaNorm builds its biology and library-size functions
from, exported so that a downstream model can use the same smooth
functions of position. It is the tensor product of a natural cubic
spline basis ([`splines::ns()`](https://rdrr.io/r/splines/ns.html))
along each axis, with its columns centred.

## Usage

``` r
tpsBasis(x, y, df = 6, ref.x = x, ref.y = y, center = TRUE)
```

## Arguments

- x, y:

  numeric vectors of the coordinates to evaluate the basis at.

- df:

  the degrees of freedom. A single positive integer applies SpaNorm's
  aspect-ratio rule to the reference ranges: with
  `gap = max(range_x, range_y) / df`, the axes get
  `ceiling(range_x / gap)` and `ceiling(range_y / gap)` degrees of
  freedom, so the longer axis gets `df`. Two positive integers give the
  x and y degrees of freedom directly. The basis has `df_x * df_y`
  columns (for example 9 at `df = c(3, 3)`).

- ref.x, ref.y:

  numeric vectors of the reference coordinates that define the knots,
  the per-axis degrees of freedom and the centring (for example the
  whole tissue section). Default: `x` and `y`.

- center:

  logical; subtract the reference column means (default `TRUE`).

## Value

a `length(x)` by `df_x * df_y` numeric matrix. Column
`(i - 1) * df_y + j` is the product of the `i`-th x-axis and the `j`-th
y-axis natural spline. Attributes: `"df.tps"`, the per-axis degrees of
freedom `c(df_x, df_y)`, and, when centred, `"scaled:center"`, the
subtracted reference means.

## Details

The basis is defined by a set of *reference* coordinates (by default the
evaluation coordinates themselves) and can be evaluated anywhere:

- the per-axis degrees of freedom come from the reference ranges (see
  `df`);

- each axis's interior knots sit at the quantiles of the reference
  coordinates along that axis and its boundary knots at their range,
  exactly as `splines::ns(ref, df = )` places them; a point outside the
  reference range is extrapolated linearly, as natural splines do;

- **centring**: when `center = TRUE`, every column has the mean of that
  column over the *reference* cells subtracted (the
  `scale(center = TRUE, scale = FALSE)` of the reference basis). The
  columns therefore average to zero over the reference cells, not over
  the evaluation cells: evaluated at a subset of the reference, each
  column is the same function of position that the whole reference sees,
  so two subsets of one section share one basis. The subtracted means
  are returned in the `"scaled:center"` attribute.

With the default reference (`ref.x = x`, `ref.y = y`) and a single `df`,
the result is identical to the internal basis SpaNorm's fit uses
(`bs.tps()`). Note that
[`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
rescales each axis to unit range before building its basis, so there the
aspect-ratio rule below gives `df` knots-worth of degrees of freedom on
both axes; to reproduce that on raw coordinates pass `df = c(df, df)`.

## Examples

``` r
set.seed(1)
x <- runif(200, 0, 2)
y <- runif(200, 0, 1)
B <- tpsBasis(x, y, df = 4)          # 4 x 2 degrees of freedom, 8 columns
dim(B)
#> [1] 200   8
# the same basis, built on the whole section and evaluated on one region
sub <- x < 0.5
Bs <- tpsBasis(x[sub], y[sub], df = 4, ref.x = x, ref.y = y)
all.equal(Bs, B[sub, ], check.attributes = FALSE)
#> [1] TRUE
```
