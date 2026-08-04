# Numerical Gradient of the Log-Density with Respect to the Response

Computes \\\partial \ell / \partial y\\ by central finite differences of
`distrib_pdf(..., log = TRUE)`. Powers the default
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
method for continuous distributions without a closed form.

## Usage

``` r
numerical_grad_y(distrib, y, theta, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  An object inheriting from class `"continuous_distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- h_rel:

  Numeric. Relative finite-difference step. Defaults to
  `.Machine$double.eps^(1/3)`.

## Value

A numeric vector of the same length as `y`.

## See also

[`numerical_hess_y`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md),
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)

## Examples

``` r
numerical_grad_y(gaussian_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> [1]  1  0 -1
```
