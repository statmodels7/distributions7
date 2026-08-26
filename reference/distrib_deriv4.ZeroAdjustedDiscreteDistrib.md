# Hurdle Fourth Derivatives

The factorization of the hurdle likelihood holds at every order, so this
is the third-order picture one step along: the components mixing `za`
with a parent parameter are exactly zero, the `za` part is the fourth
derivative of a binomial log-likelihood in one probability, and the
\\\theta\\ part is the parent's fourth derivative less that of
\\\log(1 - f_0)\\.

The truncation constant is where the work is. Its fourth derivative
comes from the moment-to-cumulant expansion over the fifteen partitions
of four indices, read on the ratios \\d^B f_0 / f_0\\ alone.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of counts.

- theta:

  A named list with the parent's parameters followed by `za`, the
  probability of a zero, in \\(0, 1)\\.

- expected:

  Logical of length 1. `TRUE` takes the expectation.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of fourth-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters and `za`, each a numeric vector of length
`length(y)`.

## See also

[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the wrapper;
[`distrib_deriv3.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroAdjustedDiscreteDistrib.md)
for the order below;
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
the partition sum the truncation constant uses.

## Examples

``` r
h <- zero_adjusted(poisson_distrib())
th <- list(mu = 3, za = 0.25)

names(distrib_deriv4(h, c(0, 2, 5), th))
#> [1] "mu_mu_mu_mu" "mu_mu_mu_za" "mu_mu_za_za" "mu_za_za_za" "za_za_za_za"

# Every mixed component vanishes, at this order as at the one below.
vapply(c("mu_mu_mu_za", "mu_mu_za_za", "mu_za_za_za"),
       function(k) distrib_deriv4(h, 2, th)[[k]], numeric(1))
#> mu_mu_mu_za mu_mu_za_za mu_za_za_za 
#>           0           0           0 
```
