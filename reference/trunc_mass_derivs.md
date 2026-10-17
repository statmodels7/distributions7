# Derivatives of the Truncation Constant via the Parent's CDF

Computes \\d^B Z = d^B F(U) - d^B F(L^-)\\ from the parent's cdf
derivatives, or returns `NULL` when that route is not available.
Dividing the result by \\Z\\ gives the truncated expectations that
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
and
[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
need, replacing one quadrature per component with two calls on the
parent.

## Usage

``` r
trunc_mass_derivs(distrib, theta, order)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- theta:

  A named list of the parent's parameters.

- order:

  The derivative order, an integer from 1 to 4. `NULL` is returned above
  4.

## Value

A named list of derivative components of \\Z\\, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
or `NULL` when the route is declined.

## The discrete correction

The lower endpoint is included in the truncated support, so the lower
tail subtracted in \\Z\\ is \\F(L)\\ minus the mass at \\L\\. Its
derivatives are reduced in the same way by the derivatives of that mass,
\\d^I F(L^-) = d^I F(L) - f(L) B_I\\, where \\B_I\\ is the complete Bell
polynomial in the parent's log-mass derivatives that
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
assembles.

## When NULL is returned

`NULL` is returned, and the caller falls back to quadrature, when the
order exceeds four, the highest order of a cdf derivative, and in two
further situations. One concerns accuracy and is tested by
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md).
The other concerns correctness and arises for a mixed parent with an
atom sitting exactly on the lower endpoint, because the derivative of
that mass is not the parent's own mass function.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, which both belong to the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
for the gate,
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
and
[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
for the callers at orders one and two,
[`trunc_deriv_k()`](https://statmodels7.github.io/distributions7/reference/trunc_deriv_k.md)
for orders three to five, and
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
for \\Z\\.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

Z <- distributions7:::trunc_constants(tn, theta)$Z
dZ <- distributions7:::trunc_mass_derivs(tn, theta, 1L)
unlist(dZ) / Z
#>          mu       sigma 
#>  0.08052166 -0.47667927 

# The quadrature route returns the same numbers, to machine precision.
unlist(distributions7:::trunc_score_mean_quad(tn, theta))
#>          mu       sigma 
#>  0.08052166 -0.47667927 

# A gamma parent has compiled cdf derivatives, so the route is available
# (the result is not NULL). Above order four the result is NULL.
tg <- truncated(gamma2_distrib(), lower = 0.5, upper = 5)
is.null(distributions7:::trunc_mass_derivs(tg, list(mu = 2, sigma2 = 1), 1L))
#> [1] FALSE
is.null(distributions7:::trunc_mass_derivs(tg, list(mu = 2, sigma2 = 1), 5L))
#> [1] TRUE
```
