# Excess Kurtosis of the Bernoulli Distribution

Closed form: \\\gamma_2 = \\1 - 6\mu(1-\mu)\\/\\\mu(1-\mu)\\\\, the
excess over the Gaussian. It reaches its minimum \\-2\\ at \\\mu =
1/2\\, which is the smallest excess kurtosis any distribution attains,
and diverges at either end of the unit interval.

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

A numeric vector of excess kurtoses, the length of `theta$mu`, at or
above \\-2\\.

## See also

[`skewness.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.BernoulliDistrib.md);
[`kurtosis.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.BinomialDistrib.md),
which is this divided by \\n\\;
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

## Examples

``` r
d <- bernoulli_distrib()

# The published form, written out.
all.equal(kurtosis(d, list(mu = 0.3)), (1 - 6 * 0.21) / 0.21)
#> [1] TRUE

# Minus two at one half, the lower bound for any distribution.
kurtosis(d, list(mu = c(0.1, 0.5, 0.9)))
#> [1]  5.111111 -2.000000  5.111111
```
