# Gaussian Fifth-Order Derivatives, Mean and Variance

Computes the six distinct fifth derivatives of the Gaussian log-density
with respect to \\\mu\\ and \\\sigma^2\\, in closed form from a compiled
kernel. Writing \\r = y - \mu\\ and \\v = \sigma^2\\, the three
components with three or more derivatives in \\\mu\\ are zero and
\$\$\ell^{(\mu\mu v v v)} = \dfrac{6}{v^4}, \qquad \ell^{(\mu v v v v)}
= \dfrac{24 r}{v^5}, \qquad \ell^{(v v v v v)} = -\dfrac{12}{v^5} +
\dfrac{60 r^2}{v^6}.\$\$ With `scale = "link"` the components are
composed with the link derivatives of both parameters by the same
routine that the families with a fifth-order kernel use.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled.

- scale:

  One of `"parameter"` (the default) or `"link"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of six numeric vectors, `mu_mu_mu_mu_mu`,
`mu_mu_mu_mu_sigma2`, `mu_mu_mu_sigma2_sigma2`,
`mu_mu_sigma2_sigma2_sigma2`, `mu_sigma2_sigma2_sigma2_sigma2` and
`sigma2_sigma2_sigma2_sigma2_sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`.

## See also

[`distrib_deriv4.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian2Distrib.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- gaussian2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma2 = 1.5)
d5 <- distrib_deriv5(d, y, th)
names(d5)
#> [1] "mu_mu_mu_mu_mu"                     "mu_mu_mu_mu_sigma2"                
#> [3] "mu_mu_mu_sigma2_sigma2"             "mu_mu_sigma2_sigma2_sigma2"        
#> [5] "mu_sigma2_sigma2_sigma2_sigma2"     "sigma2_sigma2_sigma2_sigma2_sigma2"
all.equal(d5$sigma2_sigma2_sigma2_sigma2_sigma2,
          -12 / 1.5^5 + 60 * (y - 0.4)^2 / 1.5^6)
#> [1] TRUE
```
