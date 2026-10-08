# Generalized Pareto Derivatives in the Response

Return the derivatives of the log-density in the response to order four,
and the mixed derivatives of orders one and two in the response and one
and two in \\(\sigma, \xi)\\, each a closed form from its own compiled
kernel. The derivative in the response is \\\ell^{(y)} = -(1 +
\xi)/(\sigma t)\\, with \\t = 1 + \xi y/\sigma\\, so none of these forms
divides by \\\xi\\; \\t\\ is formed as
[`distrib_deriv3.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md)
describes. Outside the support every component is `NaN`.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`.

- scale:

  `"parameter"` or `"link"`, for the mixed derivatives; the link scale
  is applied by the generic.

- ...:

  Unused.

## Value

A numeric vector for the derivatives in the response; a named list,
keyed by parameter or by parameter pair in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order, for the mixed derivatives.

## See also

[`distrib_gradient.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GPDDistrib.md)
for the derivatives in the parameters.

## Examples

``` r
d <- gpd_distrib()
th <- list(sigma = 1.5, xi = 0.3)
y <- c(0.3, 1.5, 6)
all.equal(distrib_grad_y(d, y, th), -(1 + 0.3) / (1.5 + 0.3 * y))
#> [1] TRUE
distrib_hess_y_hess(d, y, th)$xi_xi
#> [1]  0.3464090 -0.4045921 -0.1745479
```
