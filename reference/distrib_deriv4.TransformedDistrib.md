# Transformed Fourth Derivatives

Exactly the parent's fourth derivatives, read at the preimage \\x =
g^{-1}(y)\\, for the same reason as at third order: the Jacobian
\\\log\|dg^{-1}/dy\|\\ carries no \\\theta\\, so it contributes nothing
to any derivative in \\\theta\\ whatever the order.

Nothing here is approximated and no partition sum is taken. The
wrapper's work is entirely in the argument, not in the derivative.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- y:

  A numeric vector of observations on the transformed scale.

- theta:

  A named list of the **parent's** parameters; a transformation adds
  none.

- expected:

  Logical of length 1. Passed to the parent.

- ...:

  Passed to the parent's method.

## Value

A named list of fourth-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters, each a numeric vector of length `length(y)`.
A two-parameter parent gives five.

## See also

[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
for the wrapper;
[`distrib_deriv3.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.TransformedDistrib.md)
for the order below;
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
logn <- transformation(gaussian1_distrib(), exp_transform())
y <- c(1.2, 2.5)
th <- list(mu = 0.3, sigma = 1.1)

names(distrib_deriv4(logn, y, th))
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"
all.equal(distrib_deriv4(logn, y, th),
          distrib_deriv4(gaussian1_distrib(), log(y), th))
#> [1] TRUE
```
