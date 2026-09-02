# Is This a Fixed-Parameter Wrapper

Returns `TRUE` for a distribution produced by
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
in any of its three forms, which is how
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
of a
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)
COLLAPSES into one wrapper. The two describe the same law either way,
and nesting would leave the inner object's free set to be reconstructed
at every call.

## Usage

``` r
is_fixed(distrib)
```

## Arguments

- distrib:

  A `distrib` object.

## Value

`TRUE` for a `FixedContinuousDistrib`, a `FixedDiscreteDistrib` or a
`FixedMultivariateDistrib`, `FALSE` otherwise.

## See also

[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
which consults this, and
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
for the class.

## Examples

``` r
c(plain = distributions7:::is_fixed(gaussian1_distrib()),
  fixed = distributions7:::is_fixed(fixed(gaussian1_distrib(), mu = 0)),
  multivariate = distributions7:::is_fixed(
    fixed(mvgaussian1_distrib(2), mu1 = 0)))
#>        plain        fixed multivariate 
#>        FALSE         TRUE         TRUE 

# Which is why two calls collapse into one wrapper holding both values.
d <- fixed(fixed(gaussian1_distrib(), mu = 0), sigma = 1)
c(class = class(d)[1], held = paste(names(d@fixed_params), collapse = ", "))
#>                                    class 
#> "distributions7::FixedContinuousDistrib" 
#>                                     held 
#>                              "mu, sigma" 
```
