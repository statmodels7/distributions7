# A Constant Parameter Vector as One Number

Returns `x[1]` when `x` has several elements that are all equal and none
missing, and `x` unchanged otherwise. A compiled kernel that sums a
series for each observation computes it once for a scalar parameter and
copies the result, but repeats the same sum `length(x)` times for a
vector of equal values. A fit with no covariates on the parameters is
this case, and at a dispersion near \\10^{-13}\\ one sum of the negative
binomial expected information runs to millions of terms. The values
returned by the kernel are the same either way.

## Usage

``` r
constant_to_scalar(x)
```

## Arguments

- x:

  A numeric vector.

## Value

A numeric vector, of length one when `x` is constant.

## Examples

``` r
distributions7:::constant_to_scalar(rep(2, 5))
#> [1] 2
distributions7:::constant_to_scalar(c(2, 3))
#> [1] 2 3
```
