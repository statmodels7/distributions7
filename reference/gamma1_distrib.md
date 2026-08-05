# Gamma Distribution in Mean and Dispersion

Creates a gamma distribution object parametrized by its mean and a
dispersion, with \\\operatorname{Var}(Y) = \phi\mu^2\\.

## Usage

``` r
gamma1_distrib(link_mu = log_link(), link_phi = log_link())
```

## Arguments

- link_mu:

  Link function for \\\mu\\. Defaults to the log.

- link_phi:

  Link function for \\\phi\\. Defaults to the log.

## Value

An S7 object of class
[`Gamma1Distrib`](https://statmodels7.github.io/distributions7/reference/Gamma1Distrib.md).

## Details

This is the parametrization a generalized linear model uses: the
variance function is \\V(\mu) = \mu^2\\ and \\\phi\\ is the dispersion
that multiplies it, so the mean and the dispersion are orthogonal and
the score in \\\mu\\ is \\(y-\mu)/(\phi\mu^2)\\. The shape is \\1/\phi\\
and the rate \\1/(\phi\mu)\\.

It is the same law as
[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
which carries the mean and the *variance*: the two are related by
\\\sigma^2 = \phi\mu^2\\. They are separate families because the second
parameter is a different quantity in each, with its own interpretation,
standard error and interval.

Derivatives are closed form to fourth order, observed and expected.
Those in \\\phi\\ go through \\s = 1/\phi\\ and the one-variable chain
rule, so each polygamma function is evaluated once.

## See also

[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)

## Examples

``` r
d <- gamma1_distrib()
theta <- list(mu = 3, phi = 0.5)
distrib_pdf(d, c(1, 3, 5), theta)
#> [1] 0.22818539 0.18044704 0.07927554
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>      3.0      4.5 

# the same law as gamma2 with sigma2 = phi * mu^2
distrib_pdf(gamma2_distrib(), 3, list(mu = 3, sigma2 = 0.5 * 9))
#> [1] 0.180447
```
