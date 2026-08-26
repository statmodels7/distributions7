# The Univariate Trigamma Difference at One Degree of Freedom

Computes \\T_1(\nu) =
\tfrac12\left\[\psi'\\\left(\tfrac{\nu+1}{2}\right) -
\psi'\\\left(\tfrac{\nu}{2}\right)\right\] + \tfrac{1}{\nu^2}\\, the
odd-dimension base case of
[`mvt_T()`](https://statmodels7.github.io/distributions7/reference/mvt_T.md).
Above \\\nu = 100\\ it uses the expansion \\-1/\nu^3 + 1/\nu^5 -
3/\nu^7\\, and below it the direct difference.

## Usage

``` r
mvt_T1(nu)
```

## Arguments

- nu:

  The degrees of freedom, a numeric vector of any length, strictly
  positive. The branch is taken elementwise.

## Value

A numeric vector as long as `nu`, negative and decaying as \\-1/\nu^3\\.

## Details

The series comes from the same duplication formula this package's
`src/student_t.cpp` uses, halved and with the \\1/\nu^2\\ folded in.
Across the switch the two branches agree: at \\\nu = 99, 100, 101\\ the
values run \\-1.0305\times10^{-6}\\, \\-9.9990\times10^{-7}\\,
\\-9.7049\times10^{-7}\\ with no step.

Note the sign. \\T_1\\ is NEGATIVE, where \\T_p\\ for even \\p\\ is
non-negative; an odd dimension assembles from this base case plus a
non-negative sum.

## Notation

\\\nu\\ is the degrees of freedom and \\\psi'\\ the trigamma function.

## See also

[`mvt_T()`](https://statmodels7.github.io/distributions7/reference/mvt_T.md),
which calls this on an odd dimension, and
[`mvt_A1()`](https://statmodels7.github.io/distributions7/reference/mvt_A1.md)
for the digamma twin.

## Examples

``` r
# No step across the crossover at nu = 100.
signif(distributions7:::mvt_T1(c(99, 100, 101)), 8)
#> [1] -1.030505e-06 -9.999000e-07 -9.704950e-07

# Scaled by nu^3 the series approaches -1.
nu <- c(1e2, 1e3, 1e4)
distributions7:::mvt_T1(nu) * nu^3
#> [1] -0.999900 -0.999999 -1.000000
```
