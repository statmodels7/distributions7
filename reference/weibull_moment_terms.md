# Standardized Moments of a Weibull by Shape

Returns the pieces of the variance, the skewness and the excess kurtosis
of a Weibull with shape \\k\\, which depend only on the shape. With \\t
= 1/k\\, \\L(x) = \log\Gamma(1 + x)\\ and \\D_j = L(jt) - jL(t)\\, the
raw moments divided by \\\Gamma(1+t)^j\\ are \\e^{D_j}\\, so
\$\$\frac{\mathrm{Var}(Y)}{\mu^2\Gamma(1+t)^2} = e^{D_2} - 1,\qquad
\gamma_1 = \frac{e^{D_3} - 3e^{D_2} + 2}{(e^{D_2} - 1)^{3/2}},\qquad
\gamma_2 = \frac{e^{D_4} - 4e^{D_3} + 6e^{D_2} - 3}{(e^{D_2} - 1)^2} -
3.\$\$ For \\k \ge 5\\ the two numerators are formed from the series
\\L(x) = -\gamma x + \sum\_{m \ge 2} (-1)^m \zeta(m) x^m/m\\: the terms
linear in the \\D_j\\ are combined coefficient by coefficient, so that
the leading powers of \\t\\, which cancel exactly (to \\t^2\\ in the
skewness and \\t^3\\ in the kurtosis), are never formed; the rest are
\\e^{D} - 1 - D\\ by its own series. Otherwise the quotients are taken
in log scale, so that \\\Gamma(1 + 4/k)\\ is never formed where it
overflows.

## Usage

``` r
weibull_moment_terms(sigma)
```

## Arguments

- sigma:

  The shape, a positive numeric vector.

## Value

A list with `L1` (\\L(t)\\), `lV` (\\\log(e^{D_2} - 1)\\), `skew` and
`kurt`; each is a numeric vector of the length of `sigma`.

## See also

[`variance.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
[`skewness.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Weibull1Distrib.md)
and
[`kurtosis.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Weibull1Distrib.md),
which read it.

## Examples

``` r
# At a large shape the skewness tends to that of a minimum-Gumbel law.
distributions7:::weibull_moment_terms(c(1, 1e6))$skew
#> [1]  2.000000 -1.139541
```
