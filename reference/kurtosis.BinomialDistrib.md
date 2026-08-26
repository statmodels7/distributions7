# Excess Kurtosis of the Binomial Distribution

Closed form: \\\gamma_2 = \\1 - 6\mu(1-\mu)\\/\\n\mu(1-\mu)\\\\, the
excess over the Gaussian. It is the Bernoulli's divided by \\n\\, so it
vanishes twice as fast as the skewness and can take either sign: it is
negative over the middle of the unit interval, where the support is
bounded at both ends and the tails are lighter than a Gaussian's.

## Arguments

- x:

  A `BinomialDistrib`, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with one component, `mu` (strictly between 0 and 1), a
  numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, the length of `theta$mu`.

## See also

[`skewness.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BinomialDistrib.md);
[`kurtosis.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BernoulliDistrib.md),
the one-trial case;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

## Examples

``` r
d <- binomial_distrib(size = 10)

# The Bernoulli's divided by the trial count.
all.equal(kurtosis(d, list(mu = 0.3)),
          kurtosis(bernoulli_distrib(), list(mu = 0.3)) / 10)
#> [1] TRUE

# Negative over the middle, positive towards either end.
round(kurtosis(d, list(mu = c(0.05, 0.3, 0.5, 0.95))), 4)
#> [1]  1.5053 -0.1238 -0.2000  1.5053
```
