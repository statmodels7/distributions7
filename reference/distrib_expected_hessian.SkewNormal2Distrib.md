# Skew Normal Expected Information in the Centered Parametrization

Computes the expected second derivatives by carrying the parent's
expected information through the same congruence the observed Hessian
uses, \\J^\top E\[\ell''\] J\\ with \\J\\ the Jacobian of
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md).
The first-order term of the chain rule drops out under expectation, the
score having mean zero.

The matrix is **non-singular at zero skewness**, which the direct
parametrization's is not: there the score for \\\alpha\\ is exactly
proportional to the score for the location and the information loses a
rank. Measured at \\\mu = 0\\, \\\sigma = 1\\, the eigenvalues here tend
to 2, 1 and \\1/6\\ as \\\gamma_1 \to 0\\. Removing that singularity is
what the centered parametrization is for.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- y:

  A numeric vector. Its values do not enter the result, which is an
  expectation; only its length does, through recycling.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`. The skewness
  must not be exactly zero.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  One of `"bartlett"`, `"integrate"`, `"mc"` or `"opg"`, the strategy
  the parent uses for its own expectation. Defaults to `"bartlett"`, the
  variance of the score.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of six numeric vectors, in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order. Every entry is an expectation, so it does not depend on `y`.

## Cost, and where the digits run out

The parent's own expected information is the base class's quadrature, so
this method is a chain on top of a numerical quantity: measured at 100
observations it costs about 5.2 seconds against the parent's 2.2, where
a family that writes its information out answers in a median of 0.18
milliseconds.
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
therefore returns `FALSE` here, and `approx` is read.

The congruence is a difference of terms of size \\\gamma_1^{-2/3}\\, so
the limit is approached and then lost. Measured, the \\\gamma_1\\
component is 0.16666782 against \\1/6 = 0.16666667\\ at \\\gamma_1 =
10^{-8}\\, 0.1655 at \\10^{-10}\\, and **negative** at \\10^{-12}\\. A
fit does not visit those values, and a genuinely symmetric problem is
better posed in
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## Errors

Signals an error when any element of `gamma1` is exactly zero.

## See also

[`distrib_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal2Distrib.md)
for the observed curvature,
[`expected_hessian_exact.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.SkewNormal2Distrib.md)
for why this counts as approximated, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
e <- distrib_expected_hessian(d, 0, th)
names(e)
#> [1] "mu_mu"         "sigma_sigma"   "gamma1_gamma1" "mu_sigma"     
#> [5] "mu_gamma1"     "sigma_gamma1" 

# The information is positive definite, and stays so into symmetry, where
# the direct parametrization loses a rank.
info <- function(g) {
  e <- distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, gamma1 = g))
  M <- matrix(c(e$mu_mu, e$mu_sigma, e$mu_gamma1,
                e$mu_sigma, e$sigma_sigma, e$sigma_gamma1,
                e$mu_gamma1, e$sigma_gamma1, e$gamma1_gamma1), 3, 3)
  eigen(-M, only.values = TRUE)$values
}
rbind(gamma1_0.5 = info(0.5), gamma1_1e_6 = info(1e-6))
#>                 [,1]         [,2]          [,3]
#> gamma1_0.5  1.111280 9.468654e-18 -9.620483e-17
#> gamma1_1e_6 1.000002 1.056266e-12 -2.524355e-28

# Its own component tends to 1/6.
c(limit = 1 / 6,
  at_1e_6 = -distrib_expected_hessian(d, 0,
              list(mu = 0, sigma = 1, gamma1 = 1e-6))$gamma1_gamma1)
#>        limit      at_1e_6 
#> 1.666667e-01 2.124393e-06 
```
