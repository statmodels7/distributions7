# Mean of the Binomial Distribution

Closed form: \\E\[Y\] = n\mu\\, with \\n\\ the number of trials and
\\\mu\\ the success probability. The trial count is a property of the
object and not a parameter, so it is read from `x@size` and does not
appear in `theta`.

## Arguments

- x:

  A `BinomialDistrib`, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
  carrying the trial count in its `size` property.

- theta:

  A named list with one component, `mu` (the success probability,
  strictly between 0 and 1), a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, the length of `theta$mu`, in \\(0, n)\\.

## Notation

\\n\\ is the number of trials, held on the object, and \\\mu \in (0,1)\\
the success probability.

## See also

[`variance.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.BinomialDistrib.md);
[`mean.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.BernoulliDistrib.md),
the one-trial case;
[`mean.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.BetaBinom1Distrib.md),
the overdispersed extension;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

## Examples

``` r
d <- binomial_distrib(size = 10)

# Ten trials at probability 0.3.
mean(d, list(mu = 0.3))
#> [1] 3
```
