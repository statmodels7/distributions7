# Response Derivative of a Reparametrized Distribution

The parent's, unchanged. A reparametrization touches the parameters and
leaves the response alone, so \\\partial\log f/\partial y\\ is the
parent's at the mapped parameters and no chain rule enters.

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  The response, a numeric vector.

- theta:

  A named list of the new parameters, on the new parameter scale.

- ...:

  Passed to the parent's method.

## Value

A numeric vector, the length of `y` recycled against `theta`.

## Details

Only the continuous class registers this. A discrete family has no
response derivative, and the base class refuses one, so there is nothing
for a discrete reparametrization to delegate.

## See also

[`distrib_hess_y.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.ReparamContinuousDistrib.md)
for the second order;
[`distrib_gradient.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ReparamContinuousDistrib.md),
where the map does enter;
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

# Identical to the parent's at the mapped parameters.
all.equal(distrib_grad_y(d, c(0, 1, 2), list(mu = 1, sigma2 = 4)),
          distrib_grad_y(gaussian1_distrib(), c(0, 1, 2),
                         list(mu = 1, sigma = 2)))
#> [1] TRUE
```
