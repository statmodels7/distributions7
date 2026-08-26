# Generalized Gamma Third-Order Derivatives

Computes the ten distinct third derivatives of the generalized gamma
log-density in \\a\\, \\d\\ and \\p\\, **in closed form**, through
[`gengamma_components()`](https://statmodels7.github.io/distributions7/reference/gengamma_components.md).
The log-density is elementary apart from \\\log\Gamma(d/p)\\ and
\\\exp\\p\log(y/a)\\\\, and each of those is a univariate function of a
two-variable map, so the written-out composition
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
covers every component without forming a three-variable expansion.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead: the expected third derivatives have no closed form. That is the
one place on this page where `approx` and `nsim` are read.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- y:

  A numeric vector of positive observations. With `expected = TRUE` only
  its length is read.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `y`, all strictly positive. A
  component of length 1 is recycled.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the value at the data, computed numerically.
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  One of `"integrate"` (the default here), `"bartlett"`, `"mc"` or
  `"opg"`, the strategy
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  uses. Read only when `expected = TRUE`.

- nsim:

  A single positive integer, the sample size when `approx = "mc"`. Read
  only when `expected = TRUE`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of ten numeric vectors, `a_a_a` through `p_p_p`, each of
length `max(length(y), lengths(theta))`.

## Notation

\\\ell\\ is the log-density of one observation, \\a \> 0\\ the scale,
\\d \> 0\\ and \\p \> 0\\ the two shapes, and \\\Gamma\\ the gamma
function.

## See also

[`distrib_hessian.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GenGamma1Distrib.md)
for the order below,
[`distrib_deriv4.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GenGamma1Distrib.md)
for the order above,
[`gengamma_components()`](https://statmodels7.github.io/distributions7/reference/gengamma_components.md)
for the assembly, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
y <- c(0.6, 1.4, 3.1)
th <- list(a = 2, d = 1.5, p = 1.3)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#>  [1] "a_a_a" "a_a_d" "a_a_p" "a_d_d" "a_d_p" "a_p_p" "d_d_d" "d_d_p" "d_p_p"
#> [10] "p_p_p"

# A central difference of the Hessian reproduces the pure-scale component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(a = 2 + eps, d = 1.5, p = 1.3))$a_a
dn <- distrib_hessian(d, y, list(a = 2 - eps, d = 1.5, p = 1.3))$a_a
all.equal((up - dn) / (2 * eps), d3$a_a_a, tolerance = 1e-6)
#> [1] TRUE

# And the fully mixed component, which the two compositions never both
# contribute to.
up <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 + eps))$a_d
dn <- distrib_hessian(d, y, list(a = 2, d = 1.5, p = 1.3 - eps))$a_d
all.equal((up - dn) / (2 * eps), d3$a_d_p, tolerance = 1e-6)
#> [1] TRUE
```
