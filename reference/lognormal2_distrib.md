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

An S7 object of class `Lognormal2Distrib`, inheriting from
`continuous_distrib`, with `params` `c("mean", "var")` and `link_params`
the two links given here.

## Details

The parameters of
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
describe \\\log Y\\, so neither of them is a moment of \\Y\\. Here they
are, through \$\$\mu\_{\log} = \log\dfrac{m^2}{\sqrt{v + m^2}}, \qquad
\sigma^2\_{\log} = \log\left(1 + \dfrac{v}{m^2}\right).\$\$ Every
derivative of the log-density in \\(m, v)\\ to order five, the expected
information and its derivatives, the derivatives in the response and the
derivatives of the distribution function are closed forms, each order in
its own compiled kernel (see
[`distrib_gradient.Lognormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal2Distrib.md)).

## The distribution

\$\$f(y) = \frac{1}{y\sqrt{2\pi s^{2}}}\exp\\\left\\-\frac{(\log y -
\mu_l)^{2}}{2s^{2}}\right\\, \quad s^{2} = \log\\\left(1+\frac{v}{m^{2}}
\right)\\, \\ \mu_l = \log m - \frac{s^{2}}{2}\$\$ on \\y \in (0,
\infty)\\, with \\E\[Y\] = m\\ and \\\operatorname{Var}(Y) = v\\.

## See also

[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md);
[Lognormal2Distrib](https://statmodels7.github.io/distributions7/reference/Lognormal2Distrib.md)
for the class.

## Examples

``` r
d <- lognormal2_distrib()
theta <- list(mean = 3, var = 2)
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>        3        2 
```
