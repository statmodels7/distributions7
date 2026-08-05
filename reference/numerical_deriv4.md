# Numerical Fourth-Order Derivatives of the Log-Density

Computes the unique fourth-order partial derivatives of the log-density
by second central differences of
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md).
This powers the default
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
method for distributions without a closed-form implementation.

## Usage

``` r
numerical_deriv4(distrib, y, theta, h_rel = .Machine$double.eps^(1/4))
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
  `.Machine$double.eps^(1/4)`.

## Value

A named list of fourth-derivative component vectors, keyed as in
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)`(distrib@params, 4)`.

## Details

Each component \\\partial^4 \ell /
\partial\theta_i\partial\theta_j\partial\theta_k\partial\theta_l\\ (with
\\i \le j \le k \le l\\) is obtained as the second derivative of the
Hessian entry \\(i, j)\\ along \\(\theta_k, \theta_l)\\: a three-point
stencil when \\k = l\\, a four-point cross stencil otherwise.

## See also

[`numerical_deriv3`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)

## Examples

``` r
numerical_deriv4(gaussian1_distrib(), 0, list(mu = 0, sigma = 1))
#> $mu_mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_mu_sigma
#> [1] 0
#> 
#> $mu_mu_sigma_sigma
#> [1] -6
#> 
#> $mu_sigma_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma_sigma
#> [1] 6
#> 
```
