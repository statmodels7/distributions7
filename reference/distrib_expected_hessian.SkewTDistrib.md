# Skew t Expected Hessian and Its Derivatives

Returns the expected information, and through
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
its first and second derivatives in the parameters.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations. Its length sets the length of each
  returned component; the values themselves are not read.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`, each a
  numeric vector of length 1 or of the length of `y`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html).

- approx, nsim:

  Ignored.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  The thread count passed to the family's kernels.

## Value

A named list of numeric vectors of length `length(y)`: for
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
the components keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
for the two derivatives those keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
and
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

No component has an elementary form. The family is a location-scale
family, so every component equals its value at \\\mu = 0\\, \\\sigma =
1\\ times \\\sigma^{-k}\\, with \\k\\ the number of indices on \\\mu\\
or \\\sigma\\. The value there depends on \\(\alpha, \nu)\\ alone and is
one integral over \\z\\ of the observed derivatives against the density,
taken once per distinct pair by the exp-sinh rule of
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md).

The integral is exact to the rule's accuracy, and what it integrates
carries the family's own accuracy: every derivative in \\(\mu, \sigma,
\alpha)\\ is closed form, and every one involving \\\nu\\ comes from a
single stencil on an analytic quantity, as documented on
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md).
Against an adaptive quadrature the expected information agrees to
\\9\times10^{-12}\\; the second derivative agrees with a difference of
the first to \\10^{-3}\\ relative at \\\nu = 5\\, the components
carrying \\\nu\\ several times being differenced ones.

The tail is integrated to \\\|z\| = 10^{60}\\, which leaves a relative
\\10^{-60\nu}/\nu\\ and is negligible for \\\nu \ge 0.3\\.

`approx` and `nsim` are accepted for the generic's sake and ignored.

## See also

[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
for the construction,
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)
for the quantity this is the expectation of, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
th <- list(mu = 1, sigma = 2, alpha = 2, nu = 5)
e <- distrib_expected_hessian(d, 0, th)
vapply(e, function(v) v[1], numeric(1))
#>        mu_mu  sigma_sigma  alpha_alpha        nu_nu     mu_sigma     mu_alpha 
#> -0.444217491 -0.358152434 -0.061237108 -0.003362739 -0.146373753 -0.100397269 
#>        mu_nu  sigma_alpha     sigma_nu     alpha_nu 
#>  0.003551358  0.052256647  0.022980843 -0.003039699 

# One quadrature serves every observation sharing a shape: the location and
# the scale only rescale the result.
e2 <- distrib_expected_hessian(d, 0, list(mu = -3, sigma = 1, alpha = 2, nu = 5))
e$alpha_alpha / e2$alpha_alpha
#> [1] 1
```
