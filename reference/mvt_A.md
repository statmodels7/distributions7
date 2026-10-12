# The Digamma Difference of a Multivariate t, Without the Cancellation

Computes \\A_p(\nu) = \psi\\\left(\tfrac{\nu+p}{2}\right) -
\psi\\\left(\tfrac{\nu}{2}\right) - \tfrac{p}{\nu}\\, the part of the
score in \\\nu\\ that does not involve the data, as twice the first
derivative of
[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md).

## Usage

``` r
mvt_A(nu, p)
```

## Arguments

- nu:

  The degrees of freedom, a numeric vector of positive values. Not
  validated.

- p:

  The dimension, a single positive whole number.

## Value

A numeric vector as long as `nu`: positive at \\p = 1\\, exactly zero at
\\p = 2\\, and negative for \\p \ge 3\\.

## Details

As \\\nu\\ grows the multivariate t tends to the multivariate gaussian
and every derivative in \\\nu\\ vanishes, so \\A_p\\ is a difference of
terms agreeing to leading order: \\\psi(\tfrac{\nu+p}{2}) -
\psi(\tfrac{\nu}{2})\\ is close to \\p/\nu\\, the quantity from which it
is subtracted. In floating point the direct form loses its significant
digits as \\\nu\\ grows and eventually has the wrong sign. For even
\\p\\ the recurrence \\\psi(x+1) = \psi(x) + 1/x\\ gives a sum whose
terms carry one sign, \$\$A_p(\nu) = -\sum\_{j=0}^{p/2-1}
\frac{4j}{\nu(\nu+2j)},\$\$ exactly zero at \\p = 2\\; for odd \\p\\
[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md)
sums an asymptotic series above \\\nu = \max(20, 4p)\\.

## Notation

\\\nu\\ is the degrees of freedom, \\p\\ the dimension and \\\psi\\ the
digamma function.

## See also

[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md)
for the computation,
[`mvt_D()`](https://statmodels7.github.io/distributions7/reference/mvt_D.md)
for the part that carries the data, and
[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
for the consumer.

## Examples

``` r
# At p = 2 the quantity is exactly zero and the direct form is noise.
direct <- function(nu, p) digamma((nu + p) / 2) - digamma(nu / 2) - p / nu
rbind(exact = distributions7:::mvt_A(c(3, 10, 1e6), 2),
      direct = direct(c(3, 10, 1e6), 2))
#>                 [,1]          [,2]        [,3]
#> exact   0.000000e+00  0.000000e+00 0.00000e+00
#> direct -1.110223e-16 -5.551115e-17 2.79556e-16

# At p = 4 the two agree while the direct form has digits left, and part
# company at large nu.
nu <- c(3, 1e3, 1e6, 1e8, 1e9)
signif(rbind(exact = distributions7:::mvt_A(nu, 4),
             direct = direct(nu, 4)), 5)
#>            [,1]       [,2]        [,3]       [,4]        [,5]
#> exact  -0.26667 -3.992e-06 -4.0000e-12 -4.000e-16 -4.0000e-18
#> direct -0.26667 -3.992e-06 -3.9998e-12 -2.431e-16  3.3096e-16

# The exact form decays as -p(p - 2) / (2 nu^2), which is -4 / nu^2 at
# p = 4.
cbind(nu = nu, scaled = distributions7:::mvt_A(nu, 4) * nu^2)
#>         nu    scaled
#> [1,] 3e+00 -2.400000
#> [2,] 1e+03 -3.992016
#> [3,] 1e+06 -3.999992
#> [4,] 1e+08 -4.000000
#> [5,] 1e+09 -4.000000
```
