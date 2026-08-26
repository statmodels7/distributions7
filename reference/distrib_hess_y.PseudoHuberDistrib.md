# Pseudo-Huber Second Derivative in the Response

Computes \\\partial^2 \ell / \partial y^2\\ in closed form. With \\r =
y - \mu\\ and \\D = \sqrt{\nu + (r/\sigma)^2}\\, \$\$\dfrac{\partial^2
\ell}{\partial y^2} = -\dfrac{\nu}{\sigma^2 D^3}.\$\$ The family is a
location family, so this equals the pure-location entry of
[`distrib_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md),
with no sign change; two derivatives in \\\mu\\ carry two factors of
\\-1\\. It is **negative everywhere**, so the log-density is concave in
the response, and it decays like \\\|r\|^{-3}\\ in the tails.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length
`max(length(y), length(mu), length(sigma), length(nu))`, one value per
observation, every entry negative.

## See also

[`distrib_grad_y.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.PseudoHuberDistrib.md)
for the first derivative in the response,
[`distrib_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
for the second derivatives in the parameters, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)

# The closed form, written out.
r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
all.equal(distrib_hess_y(d, y, th), -2 / (1.2^2 * D^3))
#> [1] TRUE

# A location family, so this is the pure-location entry of the Hessian.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# A central difference of the first derivative reproduces it.
eps <- 1e-5
all.equal((distrib_grad_y(d, y + eps, th) -
           distrib_grad_y(d, y - eps, th)) / (2 * eps),
          distrib_hess_y(d, y, th), tolerance = 1e-6)
#> [1] TRUE

# Negative everywhere, decaying like |r|^-3 in the tails.
distrib_hess_y(d, 0.4 + c(1, 10, 100), th)
#> [1] -3.140246e-01 -2.299931e-03 -2.398964e-06
```
