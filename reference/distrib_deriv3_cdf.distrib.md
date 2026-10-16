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

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

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

Among the shipped univariate families the discrete ones reach these
methods, and their result is exact. Every shipped continuous family
registers its own methods: closed forms, or the compiled quadratures of
[`compiled_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf_deriv_k.md).
The continuous default, one product stencil on the distribution
function, serves families defined by the user and the wrappers, and
carries a relative error of order \\10^{-5}\\ at order 3 and \\10^{-4}\\
at order 4.

## See also

[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
and
[`cdf_scale_k()`](https://statmodels7.github.io/distributions7/reference/cdf_scale_k.md),
the two halves;
[`loc_scale_deriv_cdf_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_deriv_cdf_k.md),
the closed route that five families take instead.

## Examples

``` r
# A truncated gamma reaches this method and differences its cdf.
distrib_deriv3_cdf(truncated(gamma2_distrib(), upper = 6), 2,
                   list(mu = 2, sigma2 = 1))
#> $mu_mu_mu
#> [1] -0.7245243
#> 
#> $mu_mu_sigma2
#> [1] 0.7835818
#> 
#> $mu_sigma2_sigma2
#> [1] -0.5751564
#> 
#> $sigma2_sigma2_sigma2
#> [1] 0.03537098
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
