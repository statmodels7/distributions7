# Second-Order Truncated Moment of the Parent's Derivatives

Computes \\M\_{ij} = \mathbb{E}\_T\[H\_{ij} + s_i s_j\]\\, which is
\\d\_{ij} Z / Z\\, since \\d\_{ij} f / f = H\_{ij} + s_i s_j\\. It
enters the truncated Hessian as \\d\_{ij}\ell_T = H\_{ij}(y) - M\_{ij} +
m_i m_j\\, the subtracted part being the second derivative of \\\log
Z\\.

## Usage

``` r
trunc_M(distrib, theta)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- theta:

  A named list of the parent's parameters.

## Value

A named list of numeric vectors, one component per unordered pair of
parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

Where the parent has exact cdf derivatives the whole quantity comes from
[`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md)
at order two. Otherwise it is assembled from two quadratures,
[`trunc_hess_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_hess_mean.md)
and
[`trunc_score_prod_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_prod_mean.md),
which is also why
[`trunc_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_expected_hessian.md)
cannot use the cheap route: it needs the two pieces separately and
\\d\_{ij} Z\\ only gives their sum.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_hessian.md),
which subtracts it,
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
for the first order, and
[`trunc_hess_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_hess_mean.md)
and
[`trunc_score_prod_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_prod_mean.md)
for the two halves of the fallback.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

M <- distributions7:::trunc_M(tn, theta)
unlist(M)
#>       mu_mu sigma_sigma    mu_sigma 
#> -0.39723273  0.17501392 -0.09653031 

# It is the sum of the two quadratures, whichever route was taken.
EH <- distributions7:::trunc_hess_mean(tn, theta)
ES <- distributions7:::trunc_score_prod_mean(tn, theta)
max(abs(unlist(M) - (unlist(EH) + unlist(ES))))
#> [1] 1.276756e-15
```
