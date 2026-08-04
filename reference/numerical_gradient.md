# Numerical Gradient of the Log-Density

Computes the gradient of the log-density with respect to each parameter
by central finite differences of `distrib_pdf(..., log = TRUE)`. This
powers the default
[`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
method for distributions that do not implement an analytical gradient:
any `distrib` subclass that defines only `distrib_pdf` gets its score
function for free.

## Usage

``` r
numerical_gradient(distrib, y, theta, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters (each of length 1 or `length(y)`).

- h_rel:

  Numeric. Relative step size. Defaults to `.Machine$double.eps^(1/3)`
  (optimal for central differences).

## Value

A named list (one element per parameter) of gradient vectors.

## Details

Steps are scaled by `max(1, |theta|)` and automatically shrunk near the
boundaries of `distrib@params_bounds` so that the evaluation points
remain inside the parameter domain. Accuracy is roughly `eps^(2/3)`
(about 8 significant digits): sufficient for optimization, but slower
and less precise than an analytical implementation.

## See also

[`numerical_hessian`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md),
[`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)

## Examples

``` r
numerical_gradient(gaussian_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu
#> [1] -1  0  1
#> 
#> $sigma
#> [1]  5.500279e-11 -1.000000e+00  5.500279e-11
#> 
```
