# Elastic Net Log-CDF Derivatives

Closed form at every order from one to four. Each half of the
distribution function is a truncated Gaussian: with \\z = q - \mu\\, \\s
= \sqrt c\\ and \\x = a/\sqrt c\\ it is \\e^{w}\Phi(X)\\ below the
location and \\1 - e^{w}\Phi(X)\\ above it, for \\X = \pm sz - x\\ and a
weight \\w = -\log M(x) + x^2/2 + \mathrm{const}\\ written through the
Mills ratio that the family already carries. The four registrations are
made together by
[`register_phi_terms_cdf()`](https://statmodels7.github.io/distributions7/reference/register_phi_terms_cdf.md).

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` (any real value), `lambda`
  (positive) and `alpha` (strictly between 0 and 1), each a numeric
  vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of the requested order, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md):
three components for the gradient, six for the Hessian, ten at order 3
and fifteen at order 4.

## Separable partial derivatives

Both \\s\\ and \\x\\ are a function of \\\lambda\\ times a function of
\\\alpha\\, so their mixed partials are products of one-variable ones,
formed in the term function with a product rule in \\\alpha\\ for \\x\\,
and the one-variable ones are powers of affine functions through
[`dpow_affine()`](https://statmodels7.github.io/distributions7/reference/dpow_affine.md).

## The kink at the location

The location is the non-regular direction, as it is in the Laplace that
this family contains: the second derivative of \\F\\ in \\\mu\\ jumps at
\\q = \mu\\, and the third carries a point mass there. The formulas hold
on either side of the location, so a numerical check that straddles it
compares the arithmetic with a reference that is not valid there.

## Notation

\\\mu\\ is the location, \\\lambda \> 0\\ and \\\alpha \in (0,1)\\ the
two hyperparameters, \\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\
the two rates, \\M\\ the Mills ratio, \\\Phi\\ the standard normal
distribution function and \\F\\ the elastic net's.

## See also

[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
for the construction;
[`separable_deriv()`](https://statmodels7.github.io/distributions7/reference/separable_deriv.md)
and
[`dpow_affine()`](https://statmodels7.github.io/distributions7/reference/dpow_affine.md)
for the partials;
[`distrib_grad_cdf.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.InvGauss1Distrib.md),
the other family of this shape;
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

## Examples

``` r
d <- enet_distrib()
q <- c(-1, 0.5, 2)
th <- list(mu = 0, lambda = 1, alpha = 0.5)

# Against a central difference of the cdf, which shares no arithmetic.
fd <- numerical_cdf_deriv(d, q, th, order = 1)
max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) / unlist(fd) - 1))
#> [1] 6.65741e-10

# Three, six, ten and fifteen components as the order rises.
lengths(list(distrib_grad_cdf(d, q, th),
             distrib_hess_cdf(d, q, th),
             distrib_deriv3_cdf(d, q, th),
             distrib_deriv4_cdf(d, q, th)))
#> [1]  3  6 10 15
```
