# Is an object a torch tensor?

Is an object a torch tensor?

## Usage

``` r
is_torch_tensor(x)
```

## Arguments

- x:

  an object to test.

## Value

a logical.

## Examples

``` r
is_torch_tensor(matrix(1:4, 2, 2))
#> [1] FALSE
```
