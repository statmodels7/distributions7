# Gradient of a Reparametrized Distribution

The parent's score carried by the Jacobian of the map, \$\$\ell^{(a)} =
\sum_i \ell^{(i)}\\ \frac{\partial\theta_i}{\partial\psi_a}.\$\$ It is
exact whenever the parent's score is, and the map's partials are
whatever
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
supplies: a hand-written table where the family has one, and one stencil
per partial otherwise.

## Usage

``` r
reparam_gradient(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  The response, a numeric vector.

- theta:

  A named list of the new parameters, on the new parameter scale.

- scale:

  Either `"parameter"` (the default) or `"link"`. Handled by the
  generic, which applies the chain rule onto the link scale after this
  method has returned; the method itself always answers on the parameter
  scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per new parameter, each the length
of `y` recycled against `theta`.

## Notation

\\\ell\\ is the log-density, \\\theta\\ the parent's parameters and
\\\psi\\ the new ones.

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
for the identity;
[`distrib_hessian.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ReparamContinuousDistrib.md)
for the second order;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

## Examples

``` r
d <- reparametrize(
  gaussian1_distrib(),
  map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
  params = c("mu", "sigma2"),
  bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
  links = list(mu = linkfunctions7::identity_link(),
               sigma2 = linkfunctions7::log_link())
)
th <- list(mu = 1, sigma2 = 4)

# The same numbers as the family written out by hand, to 1e-12: the map's
# partials are differenced here and written out there.
a <- distrib_gradient(d, c(0, 1, 2), th)
b <- distrib_gradient(gaussian2_distrib(), c(0, 1, 2), th)
max(abs(unlist(a[names(b)]) - unlist(b)))
#> [1] 9.87016e-13
```
