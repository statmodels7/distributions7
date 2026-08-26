# Truncated Third Derivatives

Truncation adds no parameter and one \\\theta\\-dependent normalizing
constant: \\\ell_T = \ell - \log Z\\ with \\Z(\theta) = F(U) - F(L^-)\\.
The derivatives of \\\log Z\\ follow from the moment-to-cumulant
expansion over set partitions, read on the truncated expectations
\\\mathbb{E}\_T\[\partial^B f / f\]\\: \$\$d^I \log Z = \sum\_{\pi}
(-1)^{\|\pi\|-1}(\|\pi\|-1)! \prod\_{B \in \pi}
\mathbb{E}\_T\\\left\[\frac{\partial^B f}{f}\right\].\$\$

Each **distinct block** costs one quadrature, memoized across the
partition sum, which is why a truncated derivative is far dearer than
its parent's. Where the parent has genuinely closed cdf derivatives, or
is a lattice family whose cdf derivatives are an exact sum, that route
replaces the quadrature and is about three and a half times faster.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of observations inside the truncation interval.

- theta:

  A named list of the **parent's** parameters; truncation adds none, the
  endpoints being constants.

- expected:

  Logical of length 1. `TRUE` takes the expectation under the truncated
  law.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of third-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters, each a numeric vector of length `length(y)`.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the wrapper and the endpoint convention;
[`distrib_deriv4.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.TruncatedContinuousDistrib.md)
for the order above;
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
the cheaper route where the parent has one.

## Examples

``` r
tr <- truncated(gaussian1_distrib(), lower = -1, upper = 3)
th <- list(mu = 0.3, sigma = 1.1)

round(unlist(distrib_deriv3(tr, c(0, 1, 2), th)[["mu_mu_mu"]]), 6)
#> [1] -0.144835 -0.144835 -0.144835

# log Z is what separates it from the parent: the parent's third derivative
# in mu is identically zero and the truncated one is not.
c(truncated = distrib_deriv3(tr, 1, th)[["mu_mu_mu"]],
  parent = distrib_deriv3(gaussian1_distrib(), 1, th)[["mu_mu_mu"]])
#>  truncated     parent 
#> -0.1448346  0.0000000 
```
