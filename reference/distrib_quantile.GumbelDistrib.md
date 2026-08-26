# Gumbel Quantile Function

Computes the Gumbel quantile function in closed form, \$\$Q(p; \mu,
\sigma) = \mu - \sigma \log(-\log p).\$\$ The distribution function is
strictly increasing on the whole line, so this is its exact inverse and
the round trip through
[`distrib_cdf.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GumbelDistrib.md)
returns `p`. The median is \\\mu - \sigma\log\log 2\\, about \\\mu +
0.3665\sigma\\, and lies between the mode \\\mu\\ and the mean \\\mu +
\gamma\sigma\\.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`;
  the endpoints give `-Inf` and `Inf`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles, of length
`max(length(p), length(mu), length(sigma))`.

## See also

[`distrib_cdf.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.GumbelDistrib.md),
which this inverts;
[`distrib_rng.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.GumbelDistrib.md),
which uses the same inverse-transform identity;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
th <- list(mu = 0, sigma = 1)

# The closed form, written out.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_quantile(d, p, th), -log(-log(p)))
#> [1] TRUE

# Exact inverse: the round trip returns the probabilities it was given.
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# Mode, median and mean, in that order.
c(mode = 0, median = distrib_quantile(d, 0.5, th), mean = mean(d, th))
#>      mode    median      mean 
#> 0.0000000 0.3665129 0.5772157 

# A return level: the value exceeded once in 100 periods.
distrib_quantile(d, 1 - 1 / 100, th)
#> [1] 4.600149
```
