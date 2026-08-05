# Numerical Third-Order Derivatives of the Log-Density

Computes the unique third-order partial derivatives of the log-density
by central finite differences of
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md).
This powers the default
[`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
method for distributions without a closed-form implementation, and is
the reference used to validate the analytical kernels.

## Usage

``` r
numerical_deriv3(distrib, y, theta, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters (each of length 1 or `length(y)`).

- h_rel:

  Numeric. Relative finite-difference step. Defaults to
  `.Machine$double.eps^(1/3)`.

## Value

A named list of third-derivative component vectors, keyed as in
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)`(distrib@params, 3)`.

## Details

Each component \\\partial^3 \ell /
\partial\theta_i\partial\theta_j\partial\theta_k\\ (with \\i \le j \le
k\\) is obtained by differentiating the Hessian entry \\(i, j)\\ along
\\\theta_k\\. Steps are scaled by `max(1, |theta|)` and shrunk near
parameter-domain boundaries.

## See also

[`numerical_deriv4`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md),
[`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)

## Examples

``` r
numerical_deriv3(gaussian1_distrib(), 0, list(mu = 0, sigma = 1))
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_sigma
#> [1] 2
#> 
#> $mu_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma
#> [1] -2
#> 
```
