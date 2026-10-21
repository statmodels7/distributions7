# Binomial Log-CDF Gradient

Computes the gradient in closed form, and exactly, by the companion
identity to the Poisson's, \$\$\frac{\partial}{\partial p} P(X \le k) =
-n\\\mathrm{dbinom}(k, n-1, p).\$\$ The whole sum collapses to one
binomial mass at one fewer trial, so the function neither sums nor
differences.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
  carrying the trial count in its `size` property.

- q:

  A numeric vector of quantiles. Non-integer values are floored; values
  below zero give a derivative of zero.

- theta:

  A named list with one component, `mu` (the success probability,
  strictly between 0 and 1), a numeric vector of length 1 or of the
  length of `q`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector, `mu`, the length of `q` recycled
against `theta` and the `size` of the family. The generic signals an
error for a `size` whose length is neither 1 nor that common length.

## Notation

\\p \in (0,1)\\ is the success probability, \\n\\ the number of trials
held on the object, and \\k = \lfloor q \rfloor\\.

## See also

[`distrib_grad_cdf.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.BernoulliDistrib.md),
this identity at \\n = 1\\;
[`distrib_grad_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.discrete_distrib.md),
the general sum;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

## Examples

``` r
d <- binomial_distrib(size = 10)
q <- c(2, 5, 8)

# The component is minus n times a binomial mass at one fewer trial.
all.equal(distrib_grad_cdf(d, q, list(mu = 0.3), log = FALSE)$mu,
          -10 * dbinom(q, 9, 0.3))
#> [1] TRUE
```
