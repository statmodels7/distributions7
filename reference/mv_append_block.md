# Append One Block of Derived Quantities to Another

Concatenates two results in the shape of
[`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md),
stacking the `value`, `transform` and `block` vectors and row-binding
the two Jacobians. It is what joins a parametrization's own quantities,
from
[`mv_param_block()`](https://statmodels7.github.io/distributions7/reference/mv_param_block.md),
onto the standard deviations and correlations.

## Usage

``` r
mv_append_block(out, extra)
```

## Arguments

- out:

  A named list with `value`, `jacobian`, `transform` and `block`, as
  [`mv_derived()`](https://statmodels7.github.io/distributions7/reference/mv_derived.md)
  returns.

- extra:

  A second such list, or `NULL`, in which case `out` is returned
  unchanged. Its Jacobian must have the same number of columns as
  `out`'s.

## Value

A named list of the same four components, with the rows of `extra` after
those of `out`.

## See also

[`mv_param_block()`](https://statmodels7.github.io/distributions7/reference/mv_param_block.md)
for the usual second argument and
[`mv_derived.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.MvGaussianDistrib.md)
for the caller.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0.5)
der <- mv_derived(d, theta)

# A NULL second argument leaves the first alone.
identical(distributions7:::mv_append_block(der, NULL), der)
#> [1] TRUE

# Otherwise the rows stack and the Jacobian keeps its width.
both <- distributions7:::mv_append_block(der, der)
c(rows = length(both$value), cols = ncol(both$jacobian))
#> rows cols 
#>    6    5 
```
