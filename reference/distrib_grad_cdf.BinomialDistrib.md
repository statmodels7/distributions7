# Binomial Log-CDF Gradient

Closed form, and exact: the companion identity to the Poisson's,
\$\$\frac{\partial}{\partial p} P(X \le k) = -n\\\mathrm{dbinom}(k, n-1,
p).\$\$ The whole sum collapses to one binomial mass at one fewer trial,
so nothing is summed and nothing is differenced.

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
  strictly between 0 and 1), a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector, `mu`, the length of `q` recycled
against `theta`.

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

# Minus n times a binomial mass at one fewer trial.
all.equal(distrib_grad_cdf(d, q, list(mu = 0.3), log = FALSE)$mu,
          -10 * dbinom(q, 9, 0.3))
#> [1] TRUE
```
