# Second Derivative of the Log-Density with Respect to the Response

Computes the second derivative of the log-density with respect to the
random variable \\y\\, \\\partial^2 \ell / \partial y^2\\. Defined for
continuous distributions; distributions with a closed form provide it
directly, the others fall back to finite differences of the log-density.

## Usage

``` r
distrib_hess_y(distrib, y, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- y:

  A numeric vector of observations.

- theta:

  A named list (or named numeric vector) of distribution parameters.
  Each parameter must have length 1 or `length(y)`.

- ...:

  Additional arguments passed to the specific method.

## Value

A numeric vector of the same length as `y`.

## See also

[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)

## Examples

``` r
distrib_hess_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> [1] -1 -1 -1
```
