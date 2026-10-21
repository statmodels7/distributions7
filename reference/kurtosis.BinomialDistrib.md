# Excess Kurtosis of the Binomial Distribution

The excess kurtosis over the Gaussian has the closed form \\\gamma_2 =
\\1 - 6\mu(1-\mu)\\/\\n\mu(1-\mu)\\\\. It is the Bernoulli excess
kurtosis divided by \\n\\, so it decays like \\1/n\\, against
\\1/\sqrt{n}\\ for the skewness, and it can take either sign: it is
negative where \\\mu(1-\mu) \> 1/6\\, that is for \\\mu\\ between about
0.211 and 0.789, and positive outside that interval.

## Arguments

- x:

  A `BinomialDistrib`, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with one component, `mu` (strictly between 0 and 1), a
  numeric vector of length 1 or of the number of parameter settings.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(x@size))`.

## See also

[`skewness.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BinomialDistrib.md);
[`kurtosis.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BernoulliDistrib.md),
the one-trial case;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

## Examples

``` r
d <- binomial_distrib(size = 10)

# The excess kurtosis is the Bernoulli's divided by the trial count.
all.equal(kurtosis(d, list(mu = 0.3)),
          kurtosis(bernoulli_distrib(), list(mu = 0.3)) / 10)
#> [1] TRUE

# The values are negative over the middle and positive towards either end.
round(kurtosis(d, list(mu = c(0.05, 0.3, 0.5, 0.95))), 4)
#> [1]  1.5053 -0.1238 -0.2000  1.5053
```
