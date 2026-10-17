# Truncated Fourth Derivatives

The method extends the third-order construction by one order: \\\ell_T =
\ell - \log Z\\, and \\d^I \log Z\\ comes from the moment-to-cumulant
expansion over the fifteen partitions of four indices, where each block
is an expectation \\\mathbb{E}\_T\[\partial^B f / f\]\\ under the
truncated law.

The cost is set by the number of distinct blocks: a block that occurs in
several partitions is computed once and looked up. For the component
`mu_mu_sigma_sigma`, the fifteen partitions use eight distinct blocks,
and over all components of order four with two parameters the blocks of
one to four indices number fourteen.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of observations inside the truncation interval.

- theta:

  A named list of the parent's parameters; truncation adds none.

- expected:

  Logical of length 1. `TRUE` takes the expectation under the truncated
  law.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of fourth-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters, each a numeric vector of length `length(y)`.
A two-parameter parent gives five.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the wrapper;
[`distrib_deriv3.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.TruncatedContinuousDistrib.md)
for the order below;
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
the partition sum that both orders use.

## Examples

``` r
tr <- truncated(gaussian1_distrib(), lower = -1, upper = 3)
th <- list(mu = 0.3, sigma = 1.1)

names(distrib_deriv4(tr, c(0, 1, 2), th))
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"
round(unlist(distrib_deriv4(tr, 1, th)), 5)
#>             mu_mu_mu_mu          mu_mu_mu_sigma       mu_mu_sigma_sigma 
#>                 0.13955                 0.61531                -4.35474 
#>    mu_sigma_sigma_sigma sigma_sigma_sigma_sigma 
#>               -12.45072               -12.07624 
```
