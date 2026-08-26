# Truncated Fourth Derivatives, Discrete Parent

The third-order construction one step along, with the expectations again
taken by summation over the retained support: \\\ell_T = \ell - \log
Z\\, and \\d^I \log Z\\ from the moment-to-cumulant expansion over the
fifteen partitions of four indices.

Because every block is a finite sum, this order costs what the order
below costs times the ratio of distinct blocks, and nothing is
approximated anywhere in it.

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
  law.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of fourth-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters, each a numeric vector of length `length(y)`.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the wrapper;
[`distrib_deriv3.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.TruncatedDiscreteDistrib.md)
for the order below;
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
the partition sum both orders use.

## Examples

``` r
tr <- truncated(poisson_distrib(), lower = 1, upper = 8)
th <- list(mu = 3)

distrib_deriv4(tr, c(2, 4, 6), th)
#> $mu_mu_mu_mu
#> [1] -0.06746878 -0.21561693 -0.36376508
#> 

# A truncated lattice family with both endpoints finite has a finite
# support, so every expectation behind this is an exact sum.
sum(distrib_pdf(tr, 1:8, th))
#> [1] 1
```
