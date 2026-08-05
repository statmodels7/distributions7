# Numerical Hessian of the Log-Density

Computes the observed Hessian of the log-density with respect to the
parameters by central finite differences of
`distrib_pdf(..., log = TRUE)`. This powers the default
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
method for distributions that do not implement an analytical Hessian.

## Usage

``` r
numerical_hessian(distrib, y, theta, h_rel = .Machine$double.eps^(1/4))
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters (each of length 1 or `length(y)`).

- h_rel:

  Numeric. Relative step size. Defaults to `.Machine$double.eps^(1/4)`
  (optimal for second differences).

## Value

A named list of Hessian component vectors, in
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
order (diagonal elements first, then the upper-triangular mixed
derivatives).

## Details

Diagonal components use the three-point stencil \\(\ell(\theta+h) -
2\ell(\theta) + \ell(\theta-h))/h^2\\; mixed components use the
four-point cross stencil. Steps are scaled and clamped as in
[`numerical_gradient`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md).
Accuracy is roughly `sqrt(eps)`.

## See also

[`numerical_gradient`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md),
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)

## Examples

``` r
numerical_hessian(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu_mu
#> [1] -1 -1 -1
#> 
#> $sigma_sigma
#> [1] -2  1 -2
#> 
#> $mu_sigma
#> [1]  2  0 -2
#> 
```
