# Skew t Fourth Derivatives

Computes the thirty-five fourth derivatives of the log-density, with the
discipline of
[`distrib_deriv3.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md):
the generic construction serves every component whose Hessian entry is
closed form, and the ones it would nest are replaced by one stencil
each. \\(i, \nu, \nu, \nu)\\ becomes a third difference of the
closed-form score component \\i\\, and \\(\nu, \nu, \nu, \nu)\\ a fourth
difference of the log-density. Twenty of the thirty-five involve
\\\nu\\.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length matters.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`.

- expected:

  Logical of length 1. When `FALSE`, the default, the observed
  derivatives at `y` are returned.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read only when
  `expected = TRUE`.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of thirty-five numeric vectors, one per distinct
fourth-order component, from `mu_mu_mu_mu` to `nu_nu_nu_nu`.

## The step for the pure-nu component

A fourth difference amplifies rounding by \\h^{-4}\\, so
[`fd5_fourth()`](https://statmodels7.github.io/distributions7/reference/fd5_fourth.md)
is called at **ten times**
[`skewt_nu_step()`](https://statmodels7.github.io/distributions7/reference/skewt_nu_step.md)'s
step rather than at it. The choice is measured: at the family's base
step the per-observation noise is near \\10^{-2}\\ relative, and at ten
times that it is negligible while the \\O(h^2)\\ truncation, about
\\6\times10^{-4}\\, is what remains.

`nu_nu_nu_nu` is the least accurate quantity this family reports, and it
is also the smallest: measured on four observations at \\\nu = 6\\ it is
\\-4.5\times10^{-6}\\ while `sigma_sigma_sigma_sigma` is 127. A
**relative** comparison on a component that small is not informative,
which is why the package's battery scales its order-4 check by the size
of the whole array and reports \\3.5\times10^{-4}\\ here.

## Cost

This is the dearest method in the family: the twenty \\\nu\\ components
each cost four or five evaluations of an analytic quantity over the
whole vector. Measured at \\n = 20{,}000\\ it takes about sixteen
seconds, against sixty milliseconds for the score.

## See also

[`distrib_deriv3.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md)
for the order below,
[`fd5_fourth()`](https://statmodels7.github.io/distributions7/reference/fd5_fourth.md)
for the stencil and its \\h^{-4}\\ behavior, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
d4 <- distrib_deriv4(d, y, th)
c(components = length(d4), involving_nu = sum(grepl("nu", names(d4))))
#>   components involving_nu 
#>           35           20 

# The components span five orders of magnitude, so a relative comparison
# on the smallest of them says nothing about the largest.
s <- sort(vapply(d4, function(v) sum(abs(v)), 0), decreasing = TRUE)
s[c(1, 2, length(s) - 1, length(s))]
#> sigma_sigma_sigma_sigma          mu_mu_mu_sigma             nu_nu_nu_nu 
#>            1.268203e+02            7.633429e+01            4.736387e-03 
#>          alpha_nu_nu_nu 
#>            2.024177e-03 

# A closed-form-block component against a difference of the third order.
eps <- 1e-5
rbind(analytic = d4$mu_mu_alpha_alpha,
      numeric = (distrib_deriv3(d, y, list(mu = 0, sigma = 1,
                                           alpha = 3 + eps, nu = 6))$mu_mu_alpha -
                 distrib_deriv3(d, y, list(mu = 0, sigma = 1,
                                           alpha = 3 - eps, nu = 6))$mu_mu_alpha) /
                (2 * eps))
#>                [,1]      [,2]      [,3]         [,4]
#> analytic -0.2739461 0.9096393 0.7661879 -0.007254262
#> numeric  -0.2739428 0.9096170 0.7661854 -0.007254277
```
