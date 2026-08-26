# Generalized Gamma Quantile Function

Computes \\Q(u) = a\\\\Q\_\Gamma(u;\\ d/p)\\^{1/p}\\, inverting the same
representation the distribution function uses: \\(Y/a)^p\\ is Gamma with
shape \\d/p\\, so its quantile raised to the power \\1/p\\ and scaled by
\\a\\ is the generalized gamma's. The whole computation is one `qgamma`
call, and nothing is inverted by root finding.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  when `log.p = TRUE`. Note that the family's third parameter is also
  called `p`; here the argument is the probability, and the parameter is
  read from `theta`.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of the probabilities. All three must be
  strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is the survival probability, passed through to
  `qgamma`.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm.
  Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of the length of the recycled inputs,
strictly positive for a probability in \\(0, 1)\\.

## See also

[`distrib_cdf.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GenGamma1Distrib.md),
which it inverts,
[`distrib_rng.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GenGamma1Distrib.md)
for draws, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
th <- list(a = 2, d = 3, p = 1.5)

# The round trip through the distribution function.
y <- c(0.5, 1.5, 4)
all.equal(distrib_quantile(d, distrib_cdf(d, y, th), th), y)
#> [1] TRUE

# The identity written out.
u <- c(0.1, 0.5, 0.9)
all.equal(distrib_quantile(d, u, th),
          2 * qgamma(u, shape = 3 / 1.5)^(1 / 1.5))
#> [1] TRUE

# The median sits below the mean, the family being right-skewed here.
c(median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#>   median     mean 
#> 2.824562 3.009151 
```
