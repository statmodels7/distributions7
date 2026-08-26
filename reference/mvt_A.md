# The Digamma Difference of a Multivariate t, Without the Cancellation

Computes \\A_p(\nu) = \psi\\\left(\tfrac{\nu+p}{2}\right) -
\psi\\\left(\tfrac{\nu}{2}\right) - \tfrac{p}{\nu}\\, the part of the
score in \\\nu\\ that does not involve the data, in a form whose terms
carry one sign so that nothing cancels.

## Usage

``` r
mvt_A(nu, p)
```

## Arguments

- nu:

  The degrees of freedom, a numeric vector of any length, strictly
  positive. Nothing is validated: a non-positive value reaches
  [`base::digamma()`](https://rdrr.io/r/base/Special.html) and returns
  `NaN`.

- p:

  The dimension, a single positive whole number. Even and odd take
  different branches.

## Value

A numeric vector as long as `nu`, non-positive at every \\p\\ and
exactly zero at \\p = 2\\.

## Why the direct form fails

As \\\nu\\ grows the multivariate t tends to the multivariate gaussian
and every derivative in \\\nu\\ vanishes, so \\A_p\\ is a difference of
terms agreeing to leading order: \\\psi(\tfrac{\nu+p}{2}) -
\psi(\tfrac{\nu}{2})\\ is \\p/\nu\\, and \\p/\nu\\ is what it is
subtracted from. Measured at \\p = 4\\, the direct form is wrong by a
relative \\4.9\times10^{-5}\\ at \\\nu = 10^6\\, by 39 per cent at
\\10^8\\, and it CHANGES SIGN at \\10^9\\, reading
\\+3.3\times10^{-16}\\ where the exact value is \\-4.0\times10^{-18}\\.

## The repair needs no series for an even dimension

\\p\\ is an integer dimension, so the shift between the two digamma
arguments is a whole number of steps of the recurrence \\\psi(x+1) =
\psi(x) + 1/x\\. For even \\p\\ that gives a sum whose terms all carry
one sign: \$\$A_p(\nu) = -\sum\_{j=0}^{p/2-1}
\frac{4j}{\nu(\nu+2j)}.\$\$ It is exactly zero at \\p = 2\\, where the
direct form returns noise at \\10^{-16}\\.

For odd \\p\\ the shift is a half-integer, so the recurrence carries the
quantity onto the univariate
[`mvt_A1()`](https://statmodels7.github.io/distributions7/reference/mvt_A1.md),
which keeps a series of its own above a measured crossover.

## Notation

\\\nu\\ is the degrees of freedom, \\p\\ the dimension and \\\psi\\ the
digamma function.

## See also

[`mvt_A1()`](https://statmodels7.github.io/distributions7/reference/mvt_A1.md)
for the univariate case the odd branch reduces to,
[`mvt_T()`](https://statmodels7.github.io/distributions7/reference/mvt_T.md)
for the same treatment of the trigamma pair,
[`mvt_D()`](https://statmodels7.github.io/distributions7/reference/mvt_D.md)
for the third cancellation, and
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
# company at large nu, where the direct one changes sign.
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
