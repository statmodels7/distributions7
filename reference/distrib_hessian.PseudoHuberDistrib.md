# Pseudo-Huber Observed Hessian

Computes the six distinct second derivatives of the pseudo-Huber
log-density with respect to \\\mu\\, \\\sigma\\ and \\\nu\\, one value
per observation, in closed form. With \\r = y - \mu\\, \\D = \sqrt{\nu +
(r/\sigma)^2}\\, \\R_1 = K_1'(\sqrt{\nu})/K_1(\sqrt{\nu})\\ and \\R_2 =
K_1''(\sqrt{\nu})/K_1(\sqrt{\nu})\\, \$\$\dfrac{\partial^2
\ell}{\partial \mu^2} = -\dfrac{\nu}{\sigma^2 D^3}, \qquad
\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 3\sigma^2
r^2 D^{-1} + r^4 D^{-3}}{\sigma^6},\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \nu^2} = \dfrac{1}{4D^3} + \dfrac{1}{2\nu^2} +
\dfrac{1}{4}\left(\dfrac{R_1}{\nu^{3/2}} + \dfrac{R_1^2}{\nu} -
\dfrac{R_2}{\nu}\right),\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \mu \\
\partial \sigma} = \dfrac{-2\nu\sigma^2 r - r^3}{\sigma^2(\nu\sigma^2 +
r^2)^{3/2}}, \qquad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\nu} = -\dfrac{r}{2\sigma^2 D^3}, \qquad \dfrac{\partial^2
\ell}{\partial \sigma \\ \partial \nu} = -\dfrac{r^2}{2\sigma^3
D^3}.\$\$

The curvature in \\\mu\\ is negative at every observation, unlike a
Student t's, so the log-density is concave in the location however far
out the residual is. What redescends is the score, not the curvature's
sign.

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

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of six numeric vectors, `mu_mu`, `sigma_sigma`, `nu_nu`,
`mu_sigma`, `mu_nu` and `sigma_nu`, each of length
`max(length(y), length(mu), length(sigma), length(nu))`. The six name
the distinct entries of a symmetric \\3 \times 3\\ matrix per
observation.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location,
\\\sigma \> 0\\ the scale, \\\nu \> 0\\ the shape, \\r = y - \mu\\ and
\\K_1\\ the modified Bessel function of the second kind of order one.

## See also

[`distrib_gradient.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuberDistrib.md)
for the score,
[`distrib_expected_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PseudoHuberDistrib.md)
for the expectation of this quantity,
[`distrib_deriv3.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PseudoHuberDistrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "mu_mu"       "sigma_sigma" "nu_nu"       "mu_sigma"    "mu_nu"      
#> [6] "sigma_nu"   

# The location entry, written out, and its sign.
r <- y - 0.4; D <- sqrt(2 + (r / 1.2)^2)
all.equal(h$mu_mu, -2 / (1.2^2 * D^3))
#> [1] TRUE

# numDeriv on the summed log-density reproduces the summed matrix.
fn <- function(p)
  sum(distrib_pdf(d, y, list(mu = p[1], sigma = p[2], nu = p[3]), log = TRUE))
H <- numDeriv::hessian(fn, c(0.4, 1.2, 2))
rbind(numeric = c(H[1, 1], H[2, 2], H[3, 3], H[1, 2], H[1, 3], H[2, 3]),
      closed = vapply(h, sum, numeric(1)))
#>              mu_mu sigma_sigma       nu_nu  mu_sigma       mu_nu   sigma_nu
#> numeric -0.7771603   -2.531202 -0.01787204 0.1289325 -0.02080771 -0.2039011
#> closed  -0.7771603   -2.531202 -0.01787204 0.1289325 -0.02080771 -0.2039011

# The curvature in the location stays negative however far out y is.
distrib_hessian(d, 0.4 + c(1, 10, 100), th)$mu_mu
#> [1] -3.140246e-01 -2.299931e-03 -2.398964e-06
```
