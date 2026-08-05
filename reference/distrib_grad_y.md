# Gradient of the Log-Density with Respect to the Response

Computes the first derivative of the log-density with respect to the
random variable \\y\\ (as opposed to the parameters), \\\partial \ell /
\partial y\\. This is defined for continuous distributions;
distributions with a closed form provide it directly, the others fall
back to finite differences of the log-density.

## Usage

``` r
distrib_grad_y(distrib, y, theta, ...)
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

## Examples

``` r
distrib_grad_y(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> [1]  1  0 -1
```
