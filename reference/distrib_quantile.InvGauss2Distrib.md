# Inverse Gaussian Quantile Function in Mean and Shape

Computes the inverse Gaussian quantile function by calling
[`statmod::qinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
at `mean = mu` and `dispersion = 1/lambda`. There is no closed form: the
distribution function is elementary but not invertible in elementary
terms, so `qinvgauss()` inverts it numerically. The distribution
function is strictly increasing on \\(0, \infty)\\, so the round trip
through
[`distrib_cdf.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.InvGauss2Distrib.md)
returns `p`.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. Both must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles in \\\[0, \infty\]\\, of length
`max(length(p), length(mu), length(lambda))`.

## See also

[`distrib_cdf.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.InvGauss2Distrib.md),
which this inverts;
[`distrib_rng.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.InvGauss2Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()
th <- list(mu = 2, lambda = 3)

# A central 95 percent interval, strongly asymmetric about the mean of 2.
distrib_quantile(d, c(0.025, 0.5, 0.975), th)
#> [1] 0.4023134 1.5122507 6.4366188

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# The median falls below the mean.
c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#>   median     mean 
#> 1.512251 2.000000 
```
