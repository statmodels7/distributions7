# The Constant of a Multivariate t and Its Derivatives in nu

Computes the \\k\\-th derivative of \\c_p(\nu) =
\log\Gamma\\(\nu+p)/2\\ - \log\Gamma(\nu/2) - (p/2)\log\nu\\, for \\k\\
from 0 to 4, without the cancellation of the direct form: every
derivative vanishes as \\\nu\\ grows, while the terms of the direct form
do not.

## Usage

``` r
mvt_cp_d(nu, p, k)
```

## Arguments

- nu:

  The degrees of freedom, a numeric vector of positive values.

- p:

  The dimension, a single positive whole number.

- k:

  The order of the derivative, 0 to 4.

## Value

A numeric vector as long as `nu`.

## Details

For even \\p\\ the recurrence of the gamma function telescopes into
\$\$c_p(\nu) = -\tfrac{p}{2}\log 2 + \sum\_{j=1}^{p/2-1} \log\\\left(1 +
\tfrac{2j}{\nu}\right),\$\$ whose derivatives are
\\(-1)^{k-1}(k-1)!\\\nu^{-k}\sum_j \\(1 + 2j/\nu)^{-k} - 1\\\\, each
term formed with [`base::expm1()`](https://rdrr.io/r/base/Log.html) and
[`base::log1p()`](https://rdrr.io/r/base/Log.html). For odd \\p\\ and
\\\nu \ge \max(20, 4p)\\ the series of
[`mvt_nu_coef()`](https://statmodels7.github.io/distributions7/reference/mvt_nu_coef.md)
is summed, \$\$c_p^{(k)}(\nu) = (-1)^{k-1}\sum\_{n \ge 1} g_n
\frac{(n+k-1)!}{\nu^{n+k}}, \qquad c_p(\nu) = -\tfrac{p}{2}\log 2 -
\sum\_{n \ge 1} g_n \frac{(n-1)!}{\nu^n};\$\$ below that bound the
direct form in [`base::lgamma()`](https://rdrr.io/r/base/Special.html)
and [`base::psigamma()`](https://rdrr.io/r/base/Special.html) is used,
where it keeps its digits.

## See also

[`mvt_h_d()`](https://statmodels7.github.io/distributions7/reference/mvt_h_d.md)
for the part of the log-density that carries the data,
[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md),
which is twice its first derivative, and
[`mvt_nu_coef()`](https://statmodels7.github.io/distributions7/reference/mvt_nu_coef.md)
for the series.

## Examples

``` r
nu <- c(5, 50, 1e6)
# Against the direct form, which loses its digits as nu grows.
rbind(exact = distributions7:::mvt_cp_d(nu, 3, 0L),
      direct = lgamma((nu + 3) / 2) - lgamma(nu / 2) - 1.5 * log(nu))
#>              [,1]      [,2]     [,3]
#> exact  -0.9070803 -1.024918 -1.03972
#> direct -0.9070803 -1.024918 -1.03972
```
