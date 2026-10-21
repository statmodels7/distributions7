# Skew Normal Log-CDF Derivatives

The derivatives are in closed form at every order from one to four, in
the shape as well as in the location and the scale. With \\z =
(q-\mu)/\sigma\\ the distribution function is \\\Phi(z) -
2T(z,\alpha)\\, and Owen's \\T\\ has elementary partial derivatives in
both arguments, so the integral in its definition is differentiated away
at the first order and never has to be differentiated again. The
location and the scale then enter only through \\z\\, by the same chain
rule that the other location-scale families use.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive) and `alpha` (any
  sign); each is a numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of the requested order, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md):
three components for the gradient, six for the Hessian, ten at order 3
and fifteen at order 4.

## Accuracy

The derivatives are closed forms. A product stencil on the same
distribution function loses accuracy as the order grows, by one or two
digits per order.

## Relation to the Student t and the pseudo-Huber

The shape components of the skew normal are closed through the
elementary partial derivatives of Owen's T given above. The Student t,
the pseudo-Huber and the skew t have no closed derivative in a shape
parameter; their shape components are integrals of the density's own
derivatives, taken by the compiled rule.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\alpha\\ the
shape, \\z = (q-\mu)/\sigma\\, \\\Phi\\ the standard normal distribution
function and \\T\\ Owen's T.

## See also

[`sn_cdf_std_derivs()`](https://statmodels7.github.io/distributions7/reference/sn_cdf_std_derivs.md)
and
[`sn_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/sn_cdf_deriv_k.md)
for the construction;
[`distrib_grad_cdf.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.SkewTDistrib.md),
which integrates its shape derivatives by the compiled rule;
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## Examples

``` r
d <- skewnormal1_distrib()
q <- c(-1, 0.5, 2)
th <- list(mu = 0.3, sigma = 1.2, alpha = 2)

# The location component is exact, the density itself.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE

# At shape zero the family is Gaussian, and its location and scale
# derivatives agree with the Gaussian's.
g0 <- distrib_grad_cdf(d, q, list(mu = 0.3, sigma = 1.2, alpha = 0),
                       log = FALSE)
gg <- distrib_grad_cdf(gaussian1_distrib(), q, list(mu = 0.3, sigma = 1.2),
                       log = FALSE)
max(abs(g0$mu - gg$mu))
#> [1] 2.775558e-17

# The number of components is three, six, ten and fifteen as the order rises.
lengths(list(distrib_grad_cdf(d, q, th),
             distrib_hess_cdf(d, q, th),
             distrib_deriv3_cdf(d, q, th),
             distrib_deriv4_cdf(d, q, th)))
#> [1]  3  6 10 15
```
