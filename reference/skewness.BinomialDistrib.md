# Skewness of the Binomial Distribution

The skewness has the closed form \\\gamma_1 =
(1-2\mu)/\sqrt{n\mu(1-\mu)}\\. It is the Bernoulli skewness divided by
\\\sqrt n\\, so it vanishes as the trials accumulate, at the rate that
the central limit theorem gives. It is zero at \\\mu = 1/2\\ at every
trial count.

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

A numeric vector, of length `max(length(theta$mu), length(x@size))`.

## See also

[`kurtosis.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BinomialDistrib.md),
which vanishes as \\1/n\\;
[`skewness.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BernoulliDistrib.md),
the one-trial case;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

## Examples

``` r
d <- binomial_distrib(size = 10)

# The skewness is the Bernoulli's divided by the square root of the trial
# count.
all.equal(skewness(d, list(mu = 0.3)),
          skewness(bernoulli_distrib(), list(mu = 0.3)) / sqrt(10))
#> [1] TRUE

# The skewness is zero at one half, at any trial count.
skewness(d, list(mu = 0.5))
#> [1] 0
```
