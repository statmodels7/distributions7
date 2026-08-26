# Transformed Third Derivatives

Exactly the parent's third derivatives, read at the preimage \\x =
g^{-1}(y)\\. A transformation of the response contributes
\\\log\|dg^{-1}/dy\|\\ to the log-density, and that Jacobian does not
depend on \\\theta\\, so it vanishes from every derivative in \\\theta\\
and leaves the parent's untouched.

Measured on `transformation(gaussian1_distrib(), exp_transform())`, the
whole component list is `all.equal` to the Gaussian's at `log(y)`, at
this order and at the fourth.

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

  Logical of length 1. Passed to the parent, so the parent decides
  whether an expectation is closed or approximated.

- ...:

  Passed to the parent's method, `approx` and `nsim` among them.

## Value

A named list of third-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters, each a numeric vector of length `length(y)`.

## See also

[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
for the wrapper;
[`distrib_deriv4.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.TransformedDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
# A lognormal built by transforming a Gaussian.
logn <- transformation(gaussian1_distrib(), exp_transform())
y <- c(1.2, 2.5)
th <- list(mu = 0.3, sigma = 1.1)

distrib_deriv3(logn, y, th)
#> $mu_mu_mu
#> [1] 0 0
#> 
#> $mu_mu_sigma
#> [1] 1.50263 1.50263
#> 
#> $mu_sigma_sigma
#> [1] -0.4822558  2.5256092
#> 
#> $sigma_sigma_sigma
#> [1] -1.399446  1.327388
#> 

# It is the parent's answer at the preimage, component for component.
all.equal(distrib_deriv3(logn, y, th),
          distrib_deriv3(gaussian1_distrib(), log(y), th))
#> [1] TRUE
```
