# Lognormal Distribution in the Mean and Variance of Y

Creates a lognormal distribution object parametrized by the mean and the
variance of \\Y\\ itself, rather than of \\\log Y\\.

## Usage

``` r
lognormal2_distrib(link_mean = log_link(), link_var = log_link())
```

## Arguments

- link_mean:

  Link function for the mean. Defaults to the log.

- link_var:

  Link function for the variance. Defaults to the log.

## Value

A reparametrized distribution object.

## Details

The parameters of
[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
describe \\\log Y\\, so neither of them is a moment of \\Y\\. Here they
are, through \$\$\mu\_{\log} = \log\dfrac{m^2}{\sqrt{v + m^2}}, \qquad
\sigma^2\_{\log} = \log\left(1 + \dfrac{v}{m^2}\right)\$\$ which is the
parametrization a regression on the mean wants.

Built with
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
so every derivative to fourth order, observed and expected, is exact.

## See also

[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md),
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)

## Examples

``` r
d <- lognormal2_distrib()
theta <- list(mean = 3, var = 2)
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>        3        2 
```
