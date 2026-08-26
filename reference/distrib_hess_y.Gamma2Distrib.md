# Gamma Second Derivative in the Response, Mean and Variance

Computes the second derivative of the gamma log-density with respect to
the response, in closed form at the implied shape \\\alpha =
\mu^2/\sigma^2\\: \$\$\dfrac{\partial^2 \ell}{\partial y^2} =
-\dfrac{\alpha - 1}{y^2}.\$\$ The rate drops out, the log-density being
linear in \\y\\ apart from the \\(\alpha-1)\log y\\ term. The sign
follows the shape: the log-density is concave in the response when
\\\sigma^2 \< \mu^2\\, exactly flat when \\\sigma^2 = \mu^2\\, where the
family is exponential, and convex when \\\sigma^2 \> \mu^2\\.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of strictly positive observations. At `y = 0` the
  value is infinite unless the shape is exactly 1.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive;
  they enter only through the shape.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(sigma2))`,
one value per observation.

## See also

[`distrib_grad_y.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gamma2Distrib.md)
for the first derivative in the response,
[`distrib_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma2Distrib.md)
for the curvature in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)

all.equal(distrib_hess_y(d, y, th), -(9 / 2 - 1) / y^2)
#> [1] TRUE

# Only the shape survives, so a setting with the same mu^2/sigma2 agrees.
all.equal(distrib_hess_y(d, y, th),
          distrib_hess_y(d, y, list(mu = 6, sigma2 = 8)))
#> [1] TRUE

# Concave below sigma2 = mu^2, flat at it, convex above.
vapply(c(2, 9, 20),
       function(v) distrib_hess_y(d, 3, list(mu = 3, sigma2 = v)),
       numeric(1))
#> [1] -0.38888889  0.00000000  0.06111111
```
