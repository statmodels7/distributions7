# Third and Fourth Derivatives of the Log Distribution Function

`distrib_deriv3_cdf()` returns \\\partial^3 \log
F(q)/\partial\theta_i\partial\theta_j\partial\theta_k\\ and
`distrib_deriv4_cdf()` its fourth-order analogue, on either tail and on
the natural or the logarithmic scale. Together with
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
they complete the cdf derivative surface to fourth order.

## Usage

``` r
distrib_deriv3_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)

distrib_deriv4_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors, in which case one value is returned per setting.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default; `FALSE`
  flips the sign of every component.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default. Far into a tail the probability underflows to zero
  and the result is `-Inf` or `NaN`.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list of numeric vectors, keyed as
[`deriv_names(distrib@params, 3)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
for `distrib_deriv3_cdf()` and as `deriv_names(distrib@params, 4)` for
`distrib_deriv4_cdf()`. A two-parameter family has 4 third-order and 5
fourth-order components; a one-parameter family has one of each.

## The two routes

A discrete family uses the exact finite sum of
[`discrete_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md),
and a continuous one takes a single product stencil on its analytic
distribution function through
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md).
24 of the 42 univariate families register a closed form of their own; of
the 18 that do not, the discrete ones sum exactly and the continuous
ones (beta1, beta2, chisq, gamma1, gamma2, gengamma1 and the two von
Mises) difference.

## What consumes them

Truncation.
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
needs the derivatives of its normalizing constant \\Z = F(U) - F(L^-)\\,
and with only the first two orders available it pays one quadrature per
component at orders three and four; with these it pays two calls on the
parent instead.

## Accuracy against speed

Unusually, the closed route is the slower of the two at these orders: on
a Gaussian at 1000 quantiles it costs 7.0 ms against the stencil's 4.7
ms, because it runs a Faa di Bruno pass over the response derivatives.
It is preferred for accuracy alone, the stencil being
\\1.5\times10^{-5}\\ out at order 3 and \\1.3\times10^{-4}\\ at order 4.

## Notation

\\F\\ is the distribution function, \\S = 1 - F\\ the survival function
and \\\theta\\ the parameter on its own scale.

## See also

[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
for the two orders below;
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
the consumer;
[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md)
and
[`discrete_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md)
for the two routes.

## Examples

``` r
# Four third-order components for a two-parameter family.
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

# On the upper tail every sign flips.
distrib_deriv3_cdf(gaussian1_distrib(), 1, list(mu = 0, sigma = 1),
                   lower.tail = FALSE, log = FALSE)
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_sigma
#> [1] -0.4839414
#> 
#> $mu_sigma_sigma
#> [1] -0.4839414
#> 
#> $sigma_sigma_sigma
#> [1] 0
#> 

# One component for a one-parameter family, and the sum is exact here.
distrib_deriv4_cdf(poisson_distrib(), 3, list(mu = 2))
#> $mu_mu_mu_mu
#> [1] 0.04026212
#> 
```
