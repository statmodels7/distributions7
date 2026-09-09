# Every Skew t Derivative of an Order That Avoids the Degrees of Freedom

Returns, in closed form, the components of a given order whose indices
are all drawn from \\(\mu, \sigma, \alpha)\\: ten of the twenty at order
three and fifteen of the thirty-five at order four.

## Usage

``` r
skewt_msa_derivs(distrib, y, theta, order)
```

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`.

- order:

  A single integer, the order of differentiation.

## Value

A named list of numeric vectors, a subset of
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
at that order, holding every component free of \\\nu\\.

## Details

The multi-indices and the names come from
[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
and
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
the same enumeration, so a name is never recovered by splitting a
string.

## See also

[`skewt_msa_tower()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_tower.md)
for the derivation and
[`distrib_deriv3.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md)
for the method that reads this.

## Examples

``` r
d <- skewt_distrib()
m <- distributions7:::skewt_msa_derivs(
  d, c(-0.4, 1.2), list(mu = 0, sigma = 1, alpha = 0.7, nu = 8), 3L)
names(m)
#>  [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_mu_alpha"      
#>  [4] "mu_sigma_sigma"    "mu_sigma_alpha"    "mu_alpha_alpha"   
#>  [7] "sigma_sigma_sigma" "sigma_sigma_alpha" "sigma_alpha_alpha"
#> [10] "alpha_alpha_alpha"
```
