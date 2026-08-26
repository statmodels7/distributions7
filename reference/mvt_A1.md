# The Univariate Digamma Difference at One Degree of Freedom

Computes \\A_1(\nu) = \psi\\\left(\tfrac{\nu+1}{2}\right) -
\psi\\\left(\tfrac{\nu}{2}\right) - \tfrac{1}{\nu}\\, the odd-dimension
base case of
[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md).
The half-integer shift here leaves no recurrence to telescope, so above
\\\nu = 200\\ the expansion \\1/(2\nu^2) - 1/(4\nu^4) + 1/(2\nu^6)\\ is
used, and below it the direct difference, which still has digits there.

## Usage

``` r
mvt_A1(nu)
```

## Arguments

- nu:

  The degrees of freedom, a numeric vector of any length, strictly
  positive. The branch is taken elementwise.

## Value

A numeric vector as long as `nu`, positive and decaying as
\\1/(2\nu^2)\\.

## Details

The crossover of 200 is the one measured for the same expansion in this
package's `src/student_t.cpp`, and this is the one place in the package
where that series exists twice. The two copies are pinned against each
other in the tests. Across the switch the two branches agree to eleven
figures: at \\\nu = 199, 200, 201\\ the values run
\\1.2626\times10^{-5}\\, \\1.2500\times10^{-5}\\,
\\1.2376\times10^{-5}\\ with no step.

## Notation

\\\nu\\ is the degrees of freedom and \\\psi\\ the digamma function.

## See also

[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md),
which calls this on an odd dimension, and
[`mvt_T1()`](https://statmodels7.github.io/distributions7/reference/mvt_T1.md)
for the trigamma twin.

## Examples

``` r
# No step across the crossover at nu = 200.
signif(distributions7:::mvt_A1(c(199, 200, 201)), 10)
#> [1] 1.262578e-05 1.249984e-05 1.237578e-05

# Below it the branch is the direct difference and agrees exactly.
direct <- function(nu) digamma((nu + 1) / 2) - digamma(nu / 2) - 1 / nu
identical(distributions7:::mvt_A1(c(3, 50)), direct(c(3, 50)))
#> [1] TRUE

# Above it the series is the accurate one: scaled by 2 nu^2 both approach 1,
# and the direct form starts to wander.
nu <- c(1e3, 1e4, 1e5)
rbind(series = distributions7:::mvt_A1(nu) * 2 * nu^2,
      direct = direct(nu) * 2 * nu^2)
#>             [,1] [,2]    [,3]
#> series 0.9999995    1 1.00000
#> direct 0.9999995    1 1.00001
```
