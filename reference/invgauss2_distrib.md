# Inverse Gaussian Distribution in Mean and Shape

Creates an inverse gaussian distribution object in the classical
parametrization, the mean \\\mu\\ and the shape \\\lambda\\, with
\\\operatorname{Var}(Y) = \mu^3/\lambda\\.

## Usage

``` r
invgauss2_distrib(link_mu = log_link(), link_lambda = log_link())
```

## Arguments

- link_mu:

  Link function for \\\mu\\. Defaults to the log.

- link_lambda:

  Link function for \\\lambda\\. Defaults to the log.

## Value

An S7 object of class
[`InvGauss2Distrib`](https://statmodels7.github.io/distributions7/reference/InvGauss2Distrib.md).

## Details

The same law as
[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md),
which carries a dispersion \\\phi = 1/\lambda\\. The map between them is
one coordinate at a time, the mean being untouched, which is what makes
both sets of derivatives elementary.

The response enters every derivative linearly, so every expectation
needs only \\\mathbb{E}\[Y\] = \mu\\ and all four orders are closed
form. The mean and the shape are orthogonal.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)

## Examples

``` r
d <- invgauss2_distrib()
theta <- list(mu = 2, lambda = 3)
distrib_pdf(d, c(1, 2, 3), theta)
#> [1] 0.4749088 0.2443013 0.1173551
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#> 2.000000 2.666667 

# the same law as invgauss1 with phi = 1 / lambda
distrib_pdf(invgauss1_distrib(), 2, list(mu = 2, phi = 1 / 3))
#> [1] 0.2443013
```
