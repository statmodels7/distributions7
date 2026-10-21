# Location-Scale Distribution Function Derivatives on Either Tail

This is the body shared by the location-scale families'
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
and
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md).
Every table of
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md)
is the density times a factor, so on the log scale the factors are
multiplied by \\f/P = \exp(\log f - \log P)\\, from the family's
log-density and the log probability of the requested tail: \\f\\ and
\\P\\ are not formed where they underflow, and the upper tail is the
family's own.

## Usage

``` r
loc_scale_cdf_route(distrib, q, theta, order, lower.tail, log)
```

## Arguments

- distrib:

  A location-scale distribution object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameter values.

- order:

  The derivative order, 1 to 4.

- lower.tail, log:

  Logical of length 1.

## Value

A named list of numeric vectors, keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
at the order, or as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order 2.

## Examples

``` r
d <- logistic_distrib()
distributions7:::loc_scale_cdf_route(d, 40, list(mu = 0, sigma = 1), 1L,
                                     FALSE, TRUE)
#> $mu
#> [1] 1
#> 
#> $sigma
#> [1] 40
#> 
```
