# Student t Distribution-Function Derivatives, Standard Deviation

The derivatives of \\F(q)\\ in \\(\mu, \sigma, \nu)\\ to order four. The
family is location-scale in \\(\mu, \sigma)\\ at fixed \\\nu\\, so those
components are closed forms in the density and its derivatives; the
components in \\\nu\\, derivatives of an incomplete beta function in its
parameter, have no elementary form and are differenced, as for
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

## Arguments

- distrib:

  A `StudentT2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- lower.tail, log:

  As in
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md).

- ...:

  Unused.

## Value

A named list, one component per parameter or per multi-index.

## See also

[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md),
[`partial_loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md)
and
[`partial_loc_scale_deriv_cdf_k()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_deriv_cdf_k.md),
the bodies registered here.

## Examples

``` r
d <- student_t2_distrib()
distrib_grad_cdf(d, c(-2, 0.5), list(mu = 1, sigma = 2, nu = 6), log = FALSE)
#> $mu
#> [1] -0.0491520 -0.2219956
#> 
#> $sigma
#> [1] 0.07372800 0.05549891
#> 
#> $nu
#> [1] 0.002052216 0.003788756
#> 
```
