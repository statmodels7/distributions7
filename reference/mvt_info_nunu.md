# The Expected Information in nu of a Multivariate t

Computes \$\$I\_{\nu\nu} =
\tfrac{1}{4}\left\\\psi'\\\left(\tfrac{\nu}{2}\right) -
\psi'\\\left(\tfrac{\nu+p}{2}\right)\right\\ -
\frac{p(\nu+p+4)}{2\nu(\nu+p)(\nu+p+2)}\$\$ without the cancellation of
the direct form, whose terms are of order \\\nu^{-2}\\ while the result
is of order \\\nu^{-4}\\.

## Usage

``` r
mvt_info_nunu(nu, p)
```

## Arguments

- nu:

  The degrees of freedom, a numeric vector of positive values.

- p:

  The dimension, a single positive whole number.

## Value

A numeric vector as long as `nu`, the information, which is positive;
the expected second derivative is its negative.

## Details

Both parts are Laplace transforms: \\I\_{\nu\nu} = \int_0^\infty e^{-\nu
t}\phi(t)\\dt\\ with \$\$\phi(t) = t\\g(t) - \alpha + e^{-pt} -
\gamma\\e^{-(p+2)t}, \qquad \alpha = \frac{p+4}{2(p+2)}, \quad \gamma =
\frac{p}{2(p+2)},\$\$ \\g\\ as in
[`mvt_nu_coef()`](https://statmodels7.github.io/distributions7/reference/mvt_nu_coef.md).
The first three derivatives of \\\phi\\ at zero vanish, so \\I\_{\nu\nu}
= \sum\_{k \ge 3} \phi^{(k)}(0)/\nu^{k+1}\\; this series is summed for
\\\nu \ge \max(20, 4(p+2))\\ and the direct form is used below.

## See also

[`distrib_expected_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MvStudentTDistrib.md),
the consumer, and
[`mvt_nu_coef()`](https://statmodels7.github.io/distributions7/reference/mvt_nu_coef.md)
for the series.

## Examples

``` r
nu <- c(10, 1e4, 1e8)
direct <- (trigamma(nu / 2) - trigamma((nu + 3) / 2)) / 4 -
  3 * (nu + 7) / (2 * nu * (nu + 3) * (nu + 5))
rbind(exact = distributions7:::mvt_info_nunu(nu, 3), direct = direct)
#>                [,1]         [,2]          [,3]
#> exact  0.0006826819 1.348891e-15  1.350000e-31
#> direct 0.0006826819 1.348811e-15 -2.588087e-24
```
