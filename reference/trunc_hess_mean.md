# Truncated Mean of the Parent's Hessian

Computes \\\mathbb{E}\_T\[H\_{ij}\]\\ for every unordered pair of
parameters, by quadrature.
[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
adds it to
[`trunc_score_prod_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_prod_mean.md)
where the parent's cdf derivatives are not exact, that sum being
\\d\_{ij} Z / Z\\.

## Usage

``` r
trunc_hess_mean(distrib, theta)
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

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md),
its one caller, and
[`trunc_score_prod_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_prod_mean.md)
for the other half of the sum.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

EH <- distributions7:::trunc_hess_mean(tn, theta)
unlist(EH)
#>       mu_mu sigma_sigma    mu_sigma 
#>  -0.6944444  -0.1971907  -0.1342028 

# Against the same expectation taken by simulation.
set.seed(1)
ys <- distrib_rng(tn, 100000, theta)
H <- distrib_hessian(gaussian1_distrib(), ys, theta)
round(vapply(H, mean, numeric(1)), 3)
#>       mu_mu sigma_sigma    mu_sigma 
#>      -0.694      -0.202      -0.133 
round(unlist(EH), 3)
#>       mu_mu sigma_sigma    mu_sigma 
#>      -0.694      -0.197      -0.134 
```
