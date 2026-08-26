# Beta Second Derivative in the Response, the Shapes

Computes the second derivative of the beta log-density with respect to
the response, in closed form: \$\$\dfrac{\partial^2 \ell}{\partial y^2}
= -\dfrac{\alpha - 1}{y^2} - \dfrac{\beta - 1}{(1 - y)^2}.\$\$ It is
negative throughout, so the log-density is concave in the response,
whenever both shapes are at least one; a shape below one makes the
corresponding term positive and can turn the curvature positive near
that endpoint. At \\\alpha = \beta = 1\\ it is exactly zero, the density
being the uniform.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. An endpoint makes the
  value infinite unless the corresponding shape is exactly 1.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length
`max(length(y), length(alpha), length(beta))`, one value per
observation.

## See also

[`distrib_grad_y.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Beta2Distrib.md)
for the first derivative in the response,
[`distrib_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)

all.equal(distrib_hess_y(d, y, th),
          -(2 - 1) / y^2 - (5 - 1) / (1 - y)^2)
#> [1] TRUE

# Concave everywhere while both shapes exceed one.
all(distrib_hess_y(d, y, th) < 0)
#> [1] TRUE

# Exactly flat at alpha = beta = 1, where the beta is the uniform.
distrib_hess_y(d, y, list(alpha = 1, beta = 1))
#> [1] 0 0 0
```
