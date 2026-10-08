# Pseudo-Huber Derivatives in the Response, Standard-Deviation Parametrization

The response enters the log-density only through \\r = y - \mu\\, so
\\\partial_y = -\partial\_\mu\\.
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
returns \\-R(\nu)\\ r/(\sigma^2 D)\\ and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
returns \\-R(\nu)\\ \nu/(\sigma^2 D^3)\\;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
returns the \\\mu\\ row of the Hessian with its sign changed; the third
and fourth response derivatives are the pure-\\\mu\\ components of the
third and fourth parameter derivatives, with the sign \\(-1)^k\\. The
derivatives of the response derivatives in the parameters are the
location-scale ones of
[`partial_loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md)
and its siblings.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- scale:

  `"parameter"` or `"link"`, for
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md).

- ...:

  Unused.

## Value

A numeric vector for
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md);
a named list, one component per parameter, for
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md).

## See also

[`distrib_gradient.PseudoHuber2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuber2Distrib.md).

## Examples

``` r
d <- pseudohuber2_distrib()
th <- list(mu = 0.4, sigma = 1.5, nu = 2)
y <- c(-1, 0.5, 3)
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE
```
