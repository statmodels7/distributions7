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

A reparametrized distribution object.

## Details

The Stacy parametrization of
[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
carries a scale, a shape and a power, and exposes no mean at all, which
is awkward for a family a regression would put a linear predictor on.
Since \\\mathbb{E}\[Y\] = a\\\Gamma((d+1)/p)/\Gamma(d/p)\\, the map is
\$\$a = m\\\dfrac{\Gamma(d/p)}{\Gamma((d+1)/p)}.\$\$

## The distribution

\$\$f(y) = \frac{p\\y^{d-1}}{a^{d}\\\Gamma(d/p)}\\e^{-(y/a)^{p}}, \qquad
a = \mu\\\frac{\Gamma(d/p)}{\Gamma((d+1)/p)}\$\$ on \\y \in (0,
\infty)\\.

\$\$\mathbb{E}\[Y\] = \mu\$\$

## See also

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md),
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)

## Examples

``` r
d <- gengamma2_distrib()
theta <- list(mean = 5, d = 3, p = 1.5)
mean(d, theta)
#> [1] 5
```
