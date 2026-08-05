# Gaussian Distribution in Mean and Variance

Creates a Gaussian distribution object parametrised by its mean and its
**variance**.

## Usage

``` r
gaussian2_distrib(link_mu = identity_link(), link_sigma2 = log_link())
```

## Arguments

- link_mu:

  Link function for \\\mu\\. Defaults to the identity.

- link_sigma2:

  Link function for \\\sigma^2\\. Defaults to the log.

## Value

An S7 object of class
[`Gaussian2Distrib`](https://statmodels7.github.io/distributions7/reference/Gaussian2Distrib.md).

## Details

This is the same law as
[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
in different coordinates: \\\sigma^2\\ here is the square of the
\\\sigma\\ there. The two are separate families rather than one family
under a link, because a link changes the scale a parameter is *modelled*
on while leaving the parameter what it was, whereas here the parameter,
its interpretation, its standard error and its confidence interval are
all about the variance.

The numbering follows the literature where it has one: this
parametrisation is `NO2` in gamlss.

Derivatives are closed form to fourth order, observed and expected, and
the two parameters are orthogonal.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
[`gaussian3_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md)

## Examples

``` r
d <- gaussian2_distrib()
theta <- list(mu = 1, sigma2 = 4)
distrib_pdf(d, c(0, 1, 2), theta)
#> [1] 0.1760327 0.1994711 0.1760327
variance(d, theta)
#> [1] 4

# the same law as gaussian1 with sigma = 2
distrib_pdf(gaussian1_distrib(), 0.5, list(mu = 1, sigma = 2))
#> [1] 0.1933341
```
