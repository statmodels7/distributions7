# Transformed Score

Returns the parent's score evaluated at \\x = g^{-1}(y)\\. The Jacobian
does not depend on the parameters, so the score of the transformed model
IS the parent's; the change of variables moves where the score is read,
not what it is.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- y:

  A numeric vector of observations on the transformed scale.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The links are the parent's.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per parameter, in the parent's
order.

## Details

The transformation does not depend on \\\theta\\, so \\\ell_Y(\theta; y)
= \ell_X(\theta; g^{-1}(y)) + \log\lvert J(y)\rvert\\ has a second term
that is a CONSTANT in the parameters and differentiates away. Every
derivative in \\\theta\\ is therefore the parent's own, evaluated at \\x
= g^{-1}(y)\\, and nothing is recomputed. The response derivatives are
the exception and come from the base class, the Jacobian depending on
\\y\\.

## Notation

\\g\\ is the transformation, \\J\\ the Jacobian of its inverse and
\\\ell\\ a log-density.

## See also

[`distrib_hessian.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TransformedDistrib.md)
for the second order,
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
which is NOT the parent's, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)
set.seed(2)
y <- distrib_rng(d, 40, theta)

g <- distrib_gradient(d, y, theta)
vapply(g, sum, numeric(1))
#>        mu     sigma 
#>  4.879482 10.663691 

# It is literally the parent's score at the preimage.
identical(g, distrib_gradient(gaussian1_distrib(), log(y), theta))
#> [1] TRUE

# And it agrees with a numerical derivative of the TRANSFORMED
# log-likelihood, the Jacobian having differentiated away.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 1.321474e-09
```
