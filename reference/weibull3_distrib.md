# Weibull Distribution in the Mean

Creates a Weibull distribution object parametrised by its mean and its
shape.

## Usage

``` r
weibull3_distrib(link_mean = log_link(), link_sigma = log_link())
```

## Arguments

- link_mean:

  Link function for the mean. Defaults to the log.

- link_sigma:

  Link function for the shape. Defaults to the log.

## Value

A reparametrized distribution object.

## Details

The first parameter of
[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
is the scale and not the mean: the mean is \\\mu\\\Gamma(1 +
1/\sigma)\\. Inverting that gives the map used here, \$\$\mu =
\dfrac{m}{\Gamma(1 + 1/\sigma)},\$\$ so every derivative becomes a
derivative of the gamma function, which is why
[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
is not written this way.

The number follows gamlss, where the Weibull in the mean is `WEI3`.
Leaving `weibull2` unused is deliberate: it names a different
parametrisation there.

## See also

[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md),
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)

## Examples

``` r
d <- weibull3_distrib()
theta <- list(mean = 4, sigma = 1.7)
mean(d, theta)
#> [1] 4
```
