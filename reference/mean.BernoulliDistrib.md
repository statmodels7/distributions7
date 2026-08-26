# Mean of the Bernoulli Distribution

Closed form: \\E\[Y\] = \mu\\, the success probability. The response
takes the values 0 and 1, so the mean is the probability of the second
and every other moment is a function of it: a one-parameter family has
nothing else for them to depend on.

## Arguments

- x:

  A `BernoulliDistrib`, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- theta:

  A named list with one component, `mu` (the success probability,
  strictly between 0 and 1), a numeric vector of length 1 or `n`.
  Aligned and validated by name, so a value at or outside the unit
  interval throws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, the length of `theta$mu`, strictly inside
\\(0,1)\\.

## See also

[`variance.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.BernoulliDistrib.md);
[`mean.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.BinomialDistrib.md),
the sum of `size` of these;
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

## Examples

``` r
d <- bernoulli_distrib()

# The one parameter is the mean.
mean(d, list(mu = c(0.1, 0.5, 0.9)))
#> [1] 0.1 0.5 0.9
```
