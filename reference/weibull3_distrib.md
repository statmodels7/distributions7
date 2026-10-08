# Weibull Distribution in the Mean

Creates a Weibull distribution object parametrized by its mean and its
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

An S7 object of class `Weibull3Distrib`, inheriting from
`continuous_distrib`, with `params` `c("mean", "sigma")` and
`link_params` the two links given here.

## Details

The first parameter of
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
is the scale and not the mean: the mean is \\b\\\Gamma(1 + 1/\sigma)\\.
Here the scale is \$\$b = \dfrac{m}{\Gamma(1 + 1/\sigma)},\$\$ and every
derivative in \\(m, \sigma)\\ to order five, the expected information
and its derivatives, the derivatives in the response and the derivatives
of the distribution function are closed forms, each order in its own
compiled kernel (see
[`distrib_gradient.Weibull3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull3Distrib.md)).

The number follows gamlss, where the Weibull in the mean is `WEI3`.
Leaving `weibull2` unused is deliberate: it names a different
parametrization there.

## The distribution

\$\$f(y) = \frac{\sigma}{b}\left(\frac{y}{b}\right)^{\sigma-1}
e^{-(y/b)^{\sigma}}, \qquad b = \frac{m}{\Gamma(1+1/\sigma)}\$\$ on \\y
\in (0, \infty)\\, with \\E\[Y\] = m\\.

## References

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Journal of the Royal Statistical
Society, Series C* 54, 507-554.

## See also

[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md);
[Weibull3Distrib](https://statmodels7.github.io/distributions7/reference/Weibull3Distrib.md)
for the class.

## Examples

``` r
d <- weibull3_distrib()
theta <- list(mean = 4, sigma = 1.7)
mean(d, theta)
#> [1] 4
```
