# Third and Fourth Derivatives With Respect to the Response

\\\partial^{3}\ell/\partial y^{3}\\ and \\\partial^{4}\ell/\partial
y^{4}\\, completing the sequence begun by
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

A family whose response enters only as \\y - \mu\\ gets these from the
derivatives it already has in the location, since
\\\partial^{k}\ell/\partial y^{k} = (-1)^{k}\\
\partial^{k}\ell/\partial\mu^{k}\\; the others take one stencil of the
requested order on the log-density.

As with the orders below, a discrete family has no such derivative and
the generic rejects rather than returning a difference across the
lattice.

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
