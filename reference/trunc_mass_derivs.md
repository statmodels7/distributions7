# Derivatives of the Truncation Constant via the Parent's CDF

Computes \\d^B Z = d^B F(U) - d^B F(L^-)\\ from the parent's cdf
derivatives, or returns `NULL` when that route is not available.
Dividing the result by \\Z\\ gives the truncated expectations
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

  The derivative order, 1 or 2.

## Value

A named list of derivative components of \\Z\\, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md),
or `NULL` when the route is declined.

## The discrete correction

The lower endpoint is included in the truncated support, so what leaves
\\Z\\ is \\F(L)\\ minus the mass at \\L\\. Its derivatives lose that
mass's derivatives with it, \\d^I F(L^-) = d^I F(L) - f(L) B_I\\, where
\\B_I\\ is the complete Bell polynomial in the parent's log-mass
derivatives that
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
assembles.

## When it declines

`NULL` is returned, and the caller falls back to quadrature, in two
situations. The first is accuracy, and is decided by
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md).
The second is correctness: a mixed parent with an atom sitting exactly
on the lower endpoint, whose mass derivative is not the parent's own
mass function.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
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
for the two callers, and
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

# The same numbers the quadrature route returns, to machine precision.
unlist(distributions7:::trunc_score_mean_quad(tn, theta))
#>          mu       sigma 
#>  0.08052166 -0.47667927 

# A gamma parent has no closed-form cdf derivative, so the route declines
# and the caller integrates instead.
tg <- truncated(gamma2_distrib(), lower = 0.5, upper = 5)
is.null(distributions7:::trunc_mass_derivs(tg, list(mu = 2, sigma2 = 1), 1L))
#> [1] TRUE
```
