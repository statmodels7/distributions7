# Second Response Derivative of a Reparametrized Distribution

The parent's, unchanged, for the same reason as the first: a
reparametrization moves the parameters and leaves the response where it
was. Only the continuous class registers it.

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

## See also

[`distrib_grad_y.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ReparamContinuousDistrib.md)
for the first order;
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

# A Gaussian's is -1/sigma^2, whatever the parametrization.
distrib_hess_y(d, c(0, 1, 2), list(mu = 1, sigma2 = 4))
#> [1] -0.25 -0.25 -0.25
```
