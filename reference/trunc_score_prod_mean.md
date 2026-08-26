# Truncated Mean of Products of Scores

Computes \\\mathbb{E}\_T\[s_i s_j\]\\ for every unordered pair of
parameters, by quadrature. Two consumers need it:
[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
adds it to
[`trunc_hess_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_hess_mean.md)
where the cdf route is unavailable, and
[`trunc_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_expected_hessian.md)
uses it whatever the parent is, the expected Hessian being
\\-\mathrm{Cov}\_T(s_i, s_j)\\.

## Usage

``` r
trunc_score_prod_mean(distrib, theta)
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

This is the quantity that keeps a truncated expected information
expensive even for a parent with closed-form cdf derivatives. Those give
\\d\_{ij} Z\\, which is \\\mathbb{E}\_T\[H\_{ij}\] + \mathbb{E}\_T\[s_i
s_j\]\\, and no rearrangement separates the two terms; the covariance
needs the second alone.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_expected_hessian.md)
and
[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md),
its two consumers, and
[`trunc_hess_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_hess_mean.md)
for the other half of the second-order moment.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

ES <- distributions7:::trunc_score_prod_mean(tn, theta)
unlist(ES)
#>       mu_mu sigma_sigma    mu_sigma 
#>  0.29721172  0.37220463  0.03767246 

# Against the same expectation taken by simulation.
set.seed(1)
ys <- distrib_rng(tn, 100000, theta)
g <- distrib_gradient(gaussian1_distrib(), ys, theta)
round(c(mu_mu = mean(g$mu^2), sigma_sigma = mean(g$sigma^2),
        mu_sigma = mean(g$mu * g$sigma)), 3)
#>       mu_mu sigma_sigma    mu_sigma 
#>       0.299       0.371       0.038 
round(unlist(ES), 3)
#>       mu_mu sigma_sigma    mu_sigma 
#>       0.297       0.372       0.038 
```
