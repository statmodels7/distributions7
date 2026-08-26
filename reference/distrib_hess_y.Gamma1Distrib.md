# Gamma Second Derivative in the Response, Mean and Dispersion

Computes the second derivative of the gamma log-density with respect to
the response, in closed form at the implied shape \\a = 1/\phi\\:
\$\$\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{a - 1}{y^2}.\$\$ The
rate drops out, the log-density being linear in \\y\\ apart from the
\\(a-1)\log y\\ term. The sign follows the shape: the log-density is
concave in the response for \\\phi \< 1\\, exactly flat at \\\phi = 1\\,
where the family is exponential, and convex for \\\phi \> 1\\.

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- y:

  A numeric vector of strictly positive observations. At `y = 0` the
  value is infinite unless the shape is exactly 1.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. `mu` is not read, the rate having
  canceled. `phi` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(phi))`, one value per
observation.

## See also

[`distrib_grad_y.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gamma1Distrib.md)
for the first derivative in the response,
[`distrib_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma1Distrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)

all.equal(distrib_hess_y(d, y, th), -(1 / 0.5 - 1) / y^2)
#> [1] TRUE

# The mean does not enter: only the shape survives the second derivative.
identical(distrib_hess_y(d, y, th),
          distrib_hess_y(d, y, list(mu = 300, phi = 0.5)))
#> [1] TRUE

# Concave below phi = 1, flat at it, convex above.
vapply(c(0.5, 1, 2), function(p) distrib_hess_y(d, 3, list(mu = 3, phi = p)),
       numeric(1))
#> [1] -0.11111111  0.00000000  0.05555556
```
