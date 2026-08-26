# The Trigamma Difference of a Multivariate t, Without the Cancellation

Computes \\T_p(\nu) =
\tfrac12\left\[\psi'\\\left(\tfrac{\nu+p}{2}\right) -
\psi'\\\left(\tfrac{\nu}{2}\right)\right\] + \tfrac{p}{\nu^2}\\, the
part of the curvature in \\\nu\\ that does not involve the data, written
so that its terms carry one sign. It is the second-order twin of
[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md)
and cancels for the same reason: every derivative in \\\nu\\ vanishes as
the family approaches the gaussian.

## Usage

``` r
mvt_T(nu, p)
```

## Arguments

- nu:

  The degrees of freedom, a numeric vector of any length, strictly
  positive. Nothing is validated.

- p:

  The dimension, a single positive whole number. Even and odd take
  different branches.

## Value

A numeric vector as long as `nu`, non-negative at every \\p\\ and
exactly zero at \\p = 2\\.

## Details

For even \\p\\ the recurrence \\\psi'(x+1) = \psi'(x) - 1/x^2\\
telescopes into \$\$T_p(\nu) = \sum\_{j=0}^{p/2-1}
\frac{8j(\nu+j)}{\nu^2(\nu+2j)^2},\$\$ exactly zero at \\p = 2\\ where
the direct form returns noise. For odd \\p\\ the half-integer shift
carries the quantity onto the univariate
[`mvt_T1()`](https://statmodels7.github.io/distributions7/reference/mvt_T1.md).
Measured at \\p = 4\\, the direct form is wrong by a relative
\\4.4\times10^{-5}\\ at \\\nu = 10^6\\ and by a factor of five at
\\10^8\\.

## Notation

\\\nu\\ is the degrees of freedom, \\p\\ the dimension and \\\psi'\\ the
trigamma function.

## See also

[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md)
for the first-order twin,
[`mvt_T1()`](https://statmodels7.github.io/distributions7/reference/mvt_T1.md)
for the univariate base case, and
[`distrib_hessian.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvStudentTDistrib.md)
for the consumer.

## Examples

``` r
# Exactly zero at p = 2, where the direct form is noise.
direct <- function(nu, p)
  0.5 * (trigamma((nu + p) / 2) - trigamma(nu / 2)) + p / nu^2
rbind(exact = distributions7:::mvt_T(c(3, 10, 1e6), 2),
      direct = direct(c(3, 10, 1e6), 2))
#>        [,1]          [,2]          [,3]
#> exact     0  0.000000e+00  0.000000e+00
#> direct    0 -3.469447e-18 -5.520744e-23

# At p = 4 the two agree until the cancellation bites.
nu <- c(3, 1e4, 1e6, 1e8)
signif(rbind(exact = distributions7:::mvt_T(nu, 4),
             direct = direct(nu, 4)), 5)
#>           [,1]       [,2]       [,3]       [,4]
#> exact  0.14222 7.9976e-12 8.0000e-18 8.0000e-24
#> direct 0.14222 7.9976e-12 8.0003e-18 1.5269e-24
```
