# Skewness of the Bernoulli Distribution

Closed form: \\\gamma_1 = (1-2\mu)/\sqrt{\mu(1-\mu)}\\. It is zero at
\\\mu = 1/2\\, positive below and negative above, and it diverges at
either end of the unit interval, where almost all the mass sits on one
of the two values.

## Arguments

- x:

  A `BernoulliDistrib`, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- theta:

  A named list with one component, `mu` (strictly between 0 and 1), a
  numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, the length of `theta$mu`.

## See also

[`kurtosis.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BernoulliDistrib.md);
[`skewness.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BinomialDistrib.md),
which is this divided by \\\sqrt{n}\\;
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

## Examples

``` r
d <- bernoulli_distrib()

# The published form, written out.
all.equal(skewness(d, list(mu = 0.3)), (1 - 0.6) / sqrt(0.3 * 0.7))
#> [1] TRUE

# Zero at one half, and it diverges towards either end.
skewness(d, list(mu = c(0.01, 0.3, 0.5, 0.7, 0.99)))
#> [1]  9.8493706  0.8728716  0.0000000 -0.8728716 -9.8493706
```
