# A Family's Own First or Second CDF Derivatives as a Table

Calls the family's
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
(`k = 1`) or
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
(`k = 2`) on the lower tail and the natural scale, and orders the result
as
[`deriv_names(distrib@params, k)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
the keying that
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md)
use.

## Usage

``` r
own_cdf_deriv(distrib, q, theta, k)
```

## Arguments

- distrib:

  An object inheriting from `continuous_distrib`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale.

- k:

  The order, 1 or 2.

## Value

A named list of numeric vectors, the derivatives of \\F\\ of order `k`.

## See also

[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
and
[`partial_loc_scale_deriv_cdf_k()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_deriv_cdf_k.md),
the two callers.

## Examples

``` r
d <- gamma2_distrib()
distributions7:::own_cdf_deriv(d, c(1, 2, 4), list(mu = 2, sigma2 = 1), 2L)
#> $mu_mu
#> [1]  0.50287104 -0.17842507 -0.07975016
#> 
#> $mu_sigma2
#> [1] -0.10667155  0.18695974 -0.04287053
#> 
#> $sigma2_sigma2
#> [1] -0.09893496 -0.01675309  0.02551909
#> 
```
