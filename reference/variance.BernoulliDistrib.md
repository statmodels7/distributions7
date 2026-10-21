# Variance of the Bernoulli Distribution

The variance has the closed form \\\operatorname{Var}(Y) = \mu(1-\mu)\\.
It is largest at \\\mu = 1/2\\, where it is \\1/4\\, and goes to zero at
either end, because the response becomes deterministic there. The value
\\1/4\\ is the largest variance that any distribution on \\\\0,1\\\\ can
have.

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

A numeric vector of variances, the length of `theta$mu`, in \\(0,
1/4\]\\.

## See also

[`mean.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.BernoulliDistrib.md);
[`variance.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.BinomialDistrib.md),
which is `size` times this;
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

## Examples

``` r
d <- bernoulli_distrib()

# The variance is largest at one half and vanishes at either end.
variance(d, list(mu = c(0.01, 0.3, 0.5, 0.7, 0.99)))
#> [1] 0.0099 0.2100 0.2500 0.2100 0.0099
```
