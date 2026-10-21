# Third and Fourth Derivatives With Respect to the Response

`distrib_deriv3_y()` and `distrib_deriv4_y()` return
\\\partial^{3}\ell/\partial y^{3}\\ and \\\partial^{4}\ell/\partial
y^{4}\\, which complete the sequence begun by
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).

## Usage

``` r
distrib_deriv3_y(distrib, y, theta, ...)

distrib_deriv4_y(distrib, y, theta, ...)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- ...:

  Passed to methods.

## Value

A numeric vector the length of `y`.

## Details

For a family whose response enters only as \\y - \mu\\, these
derivatives follow from the derivatives that the family already has in
the location, since \\\partial^{k}\ell/\partial y^{k} = (-1)^{k}\\
\partial^{k}\ell/\partial\mu^{k}\\. The other shipped families take them
from closed forms of their own, and a family that is defined only by its
density takes one stencil of the requested order on the log-density.

As with the lower orders, a discrete family has no such derivative, and
the call signals an error instead of returning a difference between
support points.

## See also

[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
[`numerical_deriv_y()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv_y.md)

## Examples

``` r
distrib_deriv3_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> [1] 0 0 0

distrib_deriv4_y(logistic_distrib(), 0.5, list(mu = 0, sigma = 1))
#> [1] 0.1927135
```
