# Default Third and Fourth Log-CDF Derivatives

The fallback both generics take when a family registers no closed form.
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
assembles the derivatives of \\F\\ of every order up to the one wanted,
by the exact sum for a discrete family and by one product stencil for a
continuous one, and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md)
puts the top order on the requested tail and scale.

## Arguments

- distrib:

  An object inheriting from `distrib`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

For
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md),
a named list of third-derivative components keyed as
[`deriv_names(distrib@params, 3)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md);
for
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md),
the fourth-order components keyed as `deriv_names(distrib@params, 4)`.
Each vector is the length of `q` recycled against `theta`.

## Details

18 of the 42 univariate families reach these methods. For the ten
discrete ones the result is exact; for the eight continuous ones (beta1,
beta2, chisq, gamma1, gamma2, gengamma1 and the two von Mises) it
carries the stencil's error, about \\1.5\times10^{-5}\\ at order 3 and
\\1.3\times10^{-4}\\ at order 4.

## See also

[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md),
the two halves;
[`loc_scale_deriv_cdf_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_deriv_cdf_k.md),
the closed route four families take instead.

## Examples

``` r
# A gamma reaches this method and differences its cdf.
distrib_deriv3_cdf(gamma2_distrib(), 2, list(mu = 2, sigma2 = 1))
#> $mu_mu_mu
#> [1] -0.7286459
#> 
#> $mu_mu_sigma2
#> [1] 0.7727148
#> 
#> $mu_sigma2_sigma2
#> [1] -0.5826829
#> 
#> $sigma2_sigma2_sigma2
#> [1] 0.04939443
#> 

# A beta-binomial reaches it too, and its sum is exact.
distrib_deriv4_cdf(betabinom1_distrib(size = 10), 4,
                   list(mu = 0.3, sigma = 0.5))
#> $mu_mu_mu_mu
#> [1] -27.82549
#> 
#> $mu_mu_mu_sigma
#> [1] 6.013565
#> 
#> $mu_mu_sigma_sigma
#> [1] -16.70488
#> 
#> $mu_sigma_sigma_sigma
#> [1] 2.986546
#> 
#> $sigma_sigma_sigma_sigma
#> [1] 7.551011
#> 
```
