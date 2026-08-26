# Lognormal Quantile Function

Computes the lognormal quantile function \$\$Q(p; \mu, \sigma^2) =
\exp\left\\\mu + \sigma\\\Phi^{-1}(p)\right\\, \qquad \sigma =
\sqrt{\sigma^2},\$\$ by calling
[`stats::qlnorm()`](https://rdrr.io/r/stats/Lognormal.html). It is
closed form, the log transformation being monotone, so the round trip
through
[`distrib_cdf.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Lognormal1Distrib.md)
returns `p` exactly. The median is \\e^{\mu}\\, which sits **below** the
mean \\e^{\mu + \sigma^2/2}\\.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning; `p = 0` gives 0 and `p = 1` gives `Inf`.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities, which is how a quantile deep in a tail is
  requested without the probability underflowing. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles in \\\[0, \infty\]\\, of length
`max(length(p), length(mu), length(sigma2))`.

## See also

[`distrib_cdf.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Lognormal1Distrib.md),
which this inverts;
[`distrib_rng.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Lognormal1Distrib.md),
which does not use it;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
th <- list(mu = 0.5, sigma2 = 0.36)

# A central 95 percent interval, multiplicatively symmetric about exp(mu).
q <- distrib_quantile(d, c(0.025, 0.5, 0.975), th)
q
#> [1] 0.5086585 1.6487213 5.3440211
c(q[2] / q[1], q[3] / q[2])
#> [1] 3.241313 3.241313

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# The median is exp(mu) and the mean is larger.
c(median = distrib_quantile(d, 0.5, th), exp_mu = exp(0.5),
  mean = mean(d, th))
#>   median   exp_mu     mean 
#> 1.648721 1.648721 1.973878 
```
