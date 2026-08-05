# Gaussian Distribution in Mean and Precision

Creates a Gaussian distribution object parametrized by its mean and its
**precision** \\\tau = 1/\sigma^2\\.

## Usage

``` r
gaussian3_distrib(link_mu = identity_link(), link_tau = log_link())
```

## Arguments

- link_mu:

  Link function for \\\mu\\. Defaults to the identity.

- link_tau:

  Link function for \\\tau\\. Defaults to the log.

## Value

An S7 object of class
[`Gaussian3Distrib`](https://statmodels7.github.io/distributions7/reference/Gaussian3Distrib.md).

## Details

This is the same law as
[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
in different coordinates, and a separate family for the same reason
[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
is: the parameter here *is* the precision, and that is what the
estimate, the standard error and the interval describe.

It is the flattest of the three parametrizations. Every third derivative
is free of the response, so the observed and the expected ones coincide,
and the only non-zero fourth derivative is \\\ell^{(\tau\tau\tau\tau)} =
-3/\tau^4\\.

The precision is the parametrization a Bayesian conjugate analysis uses,
the gamma being conjugate for \\\tau\\ at known \\\mu\\.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
[`gaussian2_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)

## Examples

``` r
d <- gaussian3_distrib()
theta <- list(mu = 1, tau = 0.25)
distrib_pdf(d, c(0, 1, 2), theta)
#> [1] 0.1760327 0.1994711 0.1760327
variance(d, theta)
#> [1] 4

# the same law as gaussian1 with sigma = 2
distrib_pdf(gaussian1_distrib(), 0.5, list(mu = 1, sigma = 2))
#> [1] 0.1933341
```
