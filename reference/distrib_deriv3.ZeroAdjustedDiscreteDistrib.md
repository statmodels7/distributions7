# Hurdle Third Derivatives

The hurdle likelihood **factorizes** into a binary part and a positive
part, so at every order the mixed components in `za` and a parent
parameter are exactly zero, and the two halves are differentiated
separately. The \\\theta\\ part is the parent's derivative less that of
the truncation constant \\\log(1 - f_0)\\, which is the mass the parent
puts at zero and the hurdle removes.

This is the structural difference from
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
and it is one a reader can check in a line: there the mixed block is not
zero, because inflation *adds* to the parent's mass at zero and no
single zero can be attributed to a mechanism.

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

  Logical of length 1. `TRUE` takes the expectation, an exact sum over
  the support for a lattice parent.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of third-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters and `za`, each a numeric vector of length
`length(y)`.

## See also

[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the wrapper and the counting rule it enforces;
[`distrib_deriv4.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ZeroAdjustedDiscreteDistrib.md)
for the order above;
[`distrib_deriv3.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroInflatedDistrib.md)
for the model this is not.

## Examples

``` r
h <- zero_adjusted(poisson_distrib())
th <- list(mu = 3, za = 0.25)

round(unlist(distrib_deriv3(h, c(0, 2, 5), th)[["mu_mu_mu"]]), 6)
#> [1] 0.000000 0.087229 0.309451

# The mixed block is exactly zero at every observation, which the
# zero-inflated wrapper's is not.
distrib_deriv3(h, c(0, 2, 5), th)[["mu_mu_za"]]
#> [1] 0 0 0
distrib_deriv3(zero_inflated(poisson_distrib()), 0,
               list(mu = 3, zi = 0.25))[["mu_mu_zi"]]
#> [1] -0.4462848
```
