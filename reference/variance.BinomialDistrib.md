# Variance of the Binomial Distribution

Closed form: \\\operatorname{Var}(Y) = n\mu(1-\mu)\\. The trials are
independent, so the variance is the trial count times the Bernoulli
variance; a sample of counts out of \\n\\ whose variance exceeds this is
overdispersed, and a beta-binomial is the family that carries the
excess.

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

A numeric vector of variances, the length of `theta$mu`.

## See also

[`mean.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.BinomialDistrib.md);
[`variance.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.BernoulliDistrib.md),
the one-trial case;
[`variance.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.BetaBinom1Distrib.md),
which inflates this;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

## Examples

``` r
d <- binomial_distrib(size = 10)

# Ten times the Bernoulli variance.
all.equal(variance(d, list(mu = 0.3)),
          10 * variance(bernoulli_distrib(), list(mu = 0.3)))
#> [1] TRUE
```
