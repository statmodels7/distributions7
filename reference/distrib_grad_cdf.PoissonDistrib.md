# Poisson Log-CDF Gradient

Closed form, and exact: the sum defining \\F\\ telescopes, so
\$\$\frac{\partial}{\partial\mu}\sum\_{j \le k} \frac{e^{-\mu}\mu^j}{j!}
= \sum\_{j \le k} \\f(j-1) - f(j)\\ = -f(k).\$\$ The sensitivity of the
distribution function to the mean is minus the mass at the last point
retained, so one density evaluation replaces the whole sum.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- q:

  A numeric vector of quantiles. Non-integer values are floored, as they
  are by the distribution function; values below zero give a derivative
  of zero.

- theta:

  A named list with one component, `mu` (positive), a numeric vector of
  length 1 or `n`.

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

\\\mu \> 0\\ is the mean, \\f\\ the mass function, \\F\\ the
distribution function and \\k = \lfloor q \rfloor\\.

## See also

[`distrib_grad_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.discrete_distrib.md),
the general sum this replaces;
[`distrib_grad_cdf.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.NegBin2Distrib.md),
whose mean component tends to this one;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

## Examples

``` r
d <- poisson_distrib()
q <- c(0, 2, 5)

# Minus the mass at the last point retained.
all.equal(distrib_grad_cdf(d, q, list(mu = 3), log = FALSE)$mu,
          -distrib_pdf(d, q, list(mu = 3)))
#> [1] TRUE

# On the log scale, and on the upper tail.
distrib_grad_cdf(d, q, list(mu = 3))$mu
#> [1] -1.0000000 -0.5294118 -0.1100543
distrib_grad_cdf(d, q, list(mu = 3), lower.tail = FALSE)$mu
#> [1] 0.0523957 0.3884153 1.2013976
```
