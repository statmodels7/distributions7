# Beta Second Derivative in the Response, Mean and Precision

Computes the second derivative of the beta log-density with respect to
the response, in closed form at the implied shapes \\\alpha = \mu\phi\\
and \\\beta = (1-\mu)\phi\\: \$\$\dfrac{\partial^2 \ell}{\partial y^2} =
-\dfrac{\alpha - 1}{y^2} - \dfrac{\beta - 1}{(1 - y)^2}.\$\$ It is
negative throughout, so the log-density is concave in the response,
whenever both shapes are at least one; a shape below one makes the
corresponding term positive and can turn the curvature positive near
that endpoint. At \\\mu = 1/2\\ and \\\phi = 2\\ both shapes are 1 and
the curvature is exactly zero, the density being the uniform.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. An endpoint makes the
  value infinite unless the corresponding shape is exactly 1.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(phi))`,
one value per observation.

## See also

[`distrib_grad_y.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Beta1Distrib.md)
for the first derivative in the response,
[`distrib_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta1Distrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)

a <- 0.4 * 5
b <- 0.6 * 5
all.equal(distrib_hess_y(d, y, th),
          -(a - 1) / y^2 - (b - 1) / (1 - y)^2)
#> [1] TRUE

# Concave everywhere while both shapes exceed one.
all(distrib_hess_y(d, y, th) < 0)
#> [1] TRUE

# Exactly flat at mu = 1/2, phi = 2, where the beta is the uniform.
distrib_hess_y(d, y, list(mu = 0.5, phi = 2))
#> [1] 0 0 0
```
