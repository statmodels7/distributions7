# Skew t Fourth Derivatives Carrying Exactly One Degree-of-Freedom Index

Returns the ten fourth-order components with exactly one index equal to
\\\nu\\, each as one five-point difference in \\\nu\\ of the
**closed-form** third derivative beside it.

## Usage

``` r
skewt_msa_nu1(distrib, y, theta, nu, h)
```

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`.

- nu:

  The degrees of freedom, of length 1.

- h:

  The step, from
  [`skewt_nu_step()`](https://statmodels7.github.io/distributions7/reference/skewt_nu_step.md).

## Value

A named list of ten numeric vectors.

## Details

This is the rule
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)'s
mixed components already follow, read one order up: one difference,
taken of an analytic quantity. The generic construction of
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)
would instead take a mixed second difference of the Hessian, which is
licensed but loses the digits a second difference costs.

The third-order component each one differentiates is found by dropping
the \\\nu\\ from the multi-index and looking the result up among the
third-order indices, so no name is parsed.

## See also

[`skewt_msa_derivs()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_derivs.md)
for the quantity being differenced and
[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
for the stencil.

## Examples

``` r
d <- skewt_distrib()
th <- list(mu = 0, sigma = 1, alpha = 0.7, nu = 8)
n1 <- distributions7:::skewt_msa_nu1(
  d, c(-0.4, 1.2), th, 8, distributions7:::skewt_nu_step(8))
names(n1)
#>  [1] "mu_mu_mu_nu"          "mu_mu_sigma_nu"       "mu_mu_alpha_nu"      
#>  [4] "mu_sigma_sigma_nu"    "mu_sigma_alpha_nu"    "mu_alpha_alpha_nu"   
#>  [7] "sigma_sigma_sigma_nu" "sigma_sigma_alpha_nu" "sigma_alpha_alpha_nu"
#> [10] "alpha_alpha_alpha_nu"
```
