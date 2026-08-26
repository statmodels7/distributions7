# Variance of the Negative Binomial Distribution in the NB1 Parametrization

Closed form: \\\operatorname{Var}(Y) = \mu(1+\theta)\\, linear in the
mean. This is what NB1 means: the variance is a fixed multiple of the
mean, where NB2 makes it \\\mu + \mu^2/\theta\\ and so quadratic. The
two are different models of overdispersion and are not
reparametrizations of each other.

## Arguments

- x:

  A `NegBin1Distrib`, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `theta` (positive),
  each a numeric vector of length 1 or `n`. The Poisson limit is
  \\\theta \to 0\\ here, the opposite end from NB2's.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$theta))`.

## Details

The distinction follows Cameron and Trivedi's numbering, and it matters
for regression: under NB1 the dispersion relative to a Poisson is the
same at every fitted mean, and under NB2 it grows with the mean.

## Notation

\\\mu \> 0\\ is the mean and \\\theta \> 0\\ the dispersion.

## References

Cameron, A. C. and Trivedi, P. K. (2013). *Regression Analysis of Count
Data*, 2nd edition. Cambridge University Press.

## See also

[`mean.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.NegBin1Distrib.md);
[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
the quadratic alternative;
[`variance.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.PoissonDistrib.md),
the \\\theta \to 0\\ limit;
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

## Examples

``` r
d <- negbin1_distrib()

# Linear in the mean: the ratio to the Poisson variance is constant.
variance(d, list(mu = c(1, 10, 100), theta = 2)) / c(1, 10, 100)
#> [1] 3 3 3

# NB2's ratio grows with the mean instead.
variance(negbin2_distrib(), list(mu = c(1, 10, 100), theta = 2)) / c(1, 10, 100)
#> [1]  1.5  6.0 51.0
```
