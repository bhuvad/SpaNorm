# NB unit deviance

The negative binomial unit deviance \\2\\y \log(y/\mu) - (y + 1/\phi)
\log((y + 1/\phi)/(\mu + 1/\phi))\\\\, elementwise over a genes x cells
matrix, on the CPU or – when `y` or `mu` is a torch tensor – on that
tensor's device.

## Usage

``` r
nbUnitDeviance(y, mu, phi)
```

## Arguments

- y:

  counts, genes x cells (matrix or torch tensor).

- mu:

  fitted means, same shape (matrix or torch tensor).

- phi:

  NB dispersion: a scalar, or one value per gene (row).

## Value

a genes x cells matrix, or a torch tensor when the input was one.

## Examples

``` r
y <- matrix(rpois(20, 3), 4, 5)
nbUnitDeviance(y, mu = y + 0.5, phi = 0.1)
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.04029267 0.03059289 0.04029267 0.08715366 0.97580328
#> [2,] 0.02418059 0.04029267 0.16700856 0.08715366 0.05634445
#> [3,] 0.04029267 0.03059289 0.04029267 0.02418059 0.97580328
#> [4,] 0.04029267 0.05634445 0.16700856 0.01967605 0.05634445
```
