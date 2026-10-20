# The Taylor Coefficients Behind a Multivariate t's Degrees of Freedom

Returns the Taylor coefficients at \\t = 0\\ of \\g(t) = (1 -
e^{-pt})/(1 - e^{-2t})\\, from which every quantity of a multivariate t
depending on \\\nu\\ only follows by a Laplace transform. With
\\c_p(\nu) = \log\Gamma\\(\nu+p)/2\\ - \log\Gamma(\nu/2) -
(p/2)\log\nu\\, \$\$c_p'(\nu) = \int_0^\infty e^{-\nu t}\\g(t) -
p/2\\\\dt,\$\$ so that \\c_p'(\nu) = \sum\_{n \ge 1} g_n\\
n!/\nu^{n+1}\\, a series that converges for even \\p\\ and is asymptotic
for odd \\p\\.

## Usage

``` r
mvt_nu_coef(p, K = 40L)
```

## Arguments

- p:

  The dimension, a single positive whole number.

- K:

  The highest power kept. Defaults to 40.

## Value

A numeric vector of length `K + 1`, the coefficients of \\t^0, \ldots,
t^K\\; the first is \\p/2\\.

## Details

The coefficients come from the division of the two power series of
\\(1 - e^{-pt})/t\\ and \\(1 - e^{-2t})/t\\, and are kept by dimension
in a cache, so each dimension is computed once per session.

## See also

[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md)
and
[`mvt_info_nunu()`](https://statmodels7.github.io/distributions7/reference/mvt_info_nunu.md),
the two consumers.

## Examples

``` r
g <- distributions7:::mvt_nu_coef(2)
# At p = 2 the ratio is identically 1.
g[1:4]
#> [1] 1 0 0 0
distributions7:::mvt_nu_coef(3)[1:4]
#> [1]  1.5000 -0.7500  0.5000 -0.1875
```
