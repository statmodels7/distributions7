# Third and Fourth Derivatives of the Log Distribution Function

\\\partial^{3}\log
F(q)/\partial\theta_i\partial\theta_j\partial\theta_k\\ and its
fourth-order analogue, on either tail.

## Usage

``` r
distrib_deriv3_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)

distrib_deriv4_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log:

  Logical; whether derivatives of the log probability are wanted.

- ...:

  Passed to methods.

## Value

A named list of derivative component vectors, keyed as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)`(distrib@params, 3)`
or `4`.

## Details

These complete the sequence begun by
[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md).
What consumes them is truncation: with only the first two orders
available,
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
pays one quadrature per component at orders three and four, and with
these it pays two calls on the parent instead.

A discrete family uses the exact finite sum and a continuous one one
product stencil on its analytic distribution function; see
[`discrete_cdf_deriv_k`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md)
and
[`numerical_cdf_deriv_k`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md).
A family with a closed form registers its own method, as at the orders
below.

## See also

[`distrib_hess_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)

## Examples

``` r
distrib_deriv3_cdf(gaussian1_distrib(), 1, list(mu = 0, sigma = 1))
#> $mu_mu_mu
#> [1] -0.2957182
#> 
#> $mu_mu_sigma
#> [1] 0.4449093
#> 
#> $mu_sigma_sigma
#> [1] 0.6103367
#> 
#> $sigma_sigma_sigma
#> [1] 0.2005643
#> 

distrib_deriv4_cdf(poisson_distrib(), 3, list(mu = 2))
#> $mu_mu_mu_mu
#> [1] 0.04026212
#> 
```
