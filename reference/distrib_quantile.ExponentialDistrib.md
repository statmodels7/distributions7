# Exponential Quantile Function

Computes the exponential quantile function \$\$Q(p; \mu) = -\mu \log(1 -
p)\$\$ by calling
[`stats::qexp()`](https://rdrr.io/r/stats/Exponential.html) at
`rate = 1/mu`. The median is \\\mu \log 2 \approx 0.693\mu\\, below the
mean, the distribution being right skewed at every parameter value.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. `p = 1` gives `Inf`; a value outside
  the range gives `NaN` with a warning.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `p`. `mu` must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\, and the quantile is then \\-\mu\log
  p\\, which is exact deep in the tail.

- log.p:

  Logical of length 1. When `TRUE` the values in `p` are read as
  logarithms of probabilities. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles in \\\[0, \infty)\\, of length
`max(length(p), length(mu))`.

## See also

[`distrib_cdf.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ExponentialDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
th <- list(mu = 2)

# The median is mu log 2, below the mean.
c(median = distrib_quantile(d, 0.5, th), mu_log2 = 2 * log(2), mean = 2)
#>   median  mu_log2     mean 
#> 1.386294 1.386294 2.000000 

# Exact inverse: the round trip returns the probabilities it was given.
p <- c(0.025, 0.5, 0.975)
all.equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p)
#> [1] TRUE

# Asked from the upper tail the quantile is -mu log p, exact however small.
distrib_quantile(d, 1e-300, th, lower.tail = FALSE)
#> [1] 1381.551
-2 * log(1e-300)
#> [1] 1381.551
```
