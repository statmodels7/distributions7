# Truncated Third Derivatives, Discrete Parent

The same construction as for a continuous parent, \\\ell_T = \ell - \log
Z\\ with \\d^I \log Z\\ from the moment-to-cumulant expansion over set
partitions, with every expectation taken by **summation over the
retained support** in place of a quadrature.

That makes it both exact and cheap. The support of a truncated lattice
family is finite whenever both endpoints are, so each block is a finite
sum of terms the parent already computes, and the memoization across the
partition sum applies as before.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of counts inside the truncation interval. Both
  endpoints are **included**.

- theta:

  A named list of the **parent's** parameters; truncation adds none.

- expected:

  Logical of length 1. `TRUE` takes the expectation under the truncated
  law, itself an exact sum.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of third-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters, each a numeric vector of length `length(y)`.
A one-parameter parent gives one.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the wrapper and the endpoint convention;
[`distrib_deriv4.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.TruncatedDiscreteDistrib.md)
for the order above;
[`distrib_deriv3.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.TruncatedContinuousDistrib.md),
where the same blocks cost a quadrature apiece.

## Examples

``` r
tr <- truncated(poisson_distrib(), lower = 1, upper = 8)
th <- list(mu = 3)

distrib_deriv3(tr, c(2, 4, 6), th)
#> $mu_mu_mu
#> [1] 0.1027039 0.2508520 0.3990002
#> 

# The parent's third derivative in mu is 2 y / mu^3 - ... and the truncated
# one differs by the derivative of log Z, which is the whole wrapper.
c(truncated = distrib_deriv3(tr, 4, th)[["mu_mu_mu"]],
  parent = distrib_deriv3(poisson_distrib(), 4, th)[["mu_mu_mu"]])
#> truncated    parent 
#> 0.2508520 0.2962963 
```
