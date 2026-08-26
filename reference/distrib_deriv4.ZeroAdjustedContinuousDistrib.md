# Zero-Adjusted Fourth Derivatives, Continuous Parent

The same separation as at third order, one step along, and for the same
reason: a continuous parent has no mass at zero, so no normalizing
constant enters and the log-likelihood is \\\log(1 - \pi) +
\ell(y;\theta)\\ away from the atom. The \\\theta\\ components are the
parent's fourth derivatives unchanged, the pure `za` component is that
of \\\log(1-\pi)\\, and everything mixing the two is exactly zero.

No partition sum is taken here at all, so this is the cheapest of the
six wrappers at this order.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations, zero included.

- theta:

  A named list with the parent's parameters followed by `za`, the
  probability of the atom, in \\(0, 1)\\.

- expected:

  Logical of length 1. `TRUE` takes the expectation over the mixed law.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of fourth-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters and `za`, each a numeric vector of length
`length(y)`. A two-parameter parent gives fifteen.

## See also

[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the wrapper;
[`distrib_deriv3.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroAdjustedContinuousDistrib.md)
for the order below;
[`distrib_deriv4.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ZeroAdjustedDiscreteDistrib.md),
where a truncation constant does appear.

## Examples

``` r
za <- zero_adjusted(gamma2_distrib())
th <- list(mu = 2, sigma2 = 0.7, za = 0.25)

length(distrib_deriv4(za, c(0, 1.5, 3), th))
#> [1] 15
all.equal(distrib_deriv4(za, c(1.5, 3), th)[["mu_mu_mu_mu"]],
          distrib_deriv4(gamma2_distrib(), c(1.5, 3),
                         list(mu = 2, sigma2 = 0.7))[["mu_mu_mu_mu"]])
#> [1] TRUE
```
