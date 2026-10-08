# Generalized Gamma Distribution in the Mean

Creates a generalized gamma distribution object whose first parameter is
the mean.

## Usage

``` r
gengamma2_distrib(
  link_mean = log_link(),
  link_d = log_link(),
  link_p = log_link()
)
```

## Arguments

- link_mean:

  Link function for the mean. Defaults to the log.

- link_d:

  Link function for the shape. Defaults to the log.

- link_p:

  Link function for the power. Defaults to the log.

## Value

An S7 object of class `GenGamma2Distrib`, inheriting from
`continuous_distrib`, with `params` `c("mean", "d", "p")` and
`link_params` the three links given here.

## Details

The Stacy parametrization of
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
carries a scale, a shape and a power, and exposes no mean at all, which
is awkward for a family a regression would put a linear predictor on.
Since \\\mathbb{E}\[Y\] = a\\\Gamma((d+1)/p)/\Gamma(d/p)\\, the scale
here is \$\$a = m\\\dfrac{\Gamma(d/p)}{\Gamma((d+1)/p)},\$\$ and every
derivative in \\(m, d, p)\\ to order five, the expected information and
its derivatives, and the derivatives in the response are closed forms,
each order in its own compiled kernel (see
[`distrib_gradient.GenGamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma2Distrib.md)).

## The distribution

\$\$f(y) = \frac{p\\y^{d-1}}{a^{d}\\\Gamma(d/p)}\\e^{-(y/a)^{p}}, \qquad
a = m\\\frac{\Gamma(d/p)}{\Gamma((d+1)/p)}\$\$ on \\y \in (0, \infty)\\,
with \\E\[Y\] = m\\.

## See also

[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md);
[GenGamma2Distrib](https://statmodels7.github.io/distributions7/reference/GenGamma2Distrib.md)
for the class.

## Examples

``` r
d <- gengamma2_distrib()
theta <- list(mean = 5, d = 3, p = 1.5)
mean(d, theta)
#> [1] 5
```
