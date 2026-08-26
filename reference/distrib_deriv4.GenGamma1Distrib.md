# Generalized Gamma Fourth-Order Derivatives

Computes the fifteen distinct fourth derivatives of the generalized
gamma log-density in \\a\\, \\d\\ and \\p\\, **in closed form**, by the
construction
[`distrib_deriv3.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
describes carried one order further: the two two-variable compositions
of
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
plus the elementary pieces of the log-density.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead: the expected fourth derivatives have no closed form. That is
the one place on this page where `approx` and `nsim` are read.

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
  `"opg"`. Read only when `expected = TRUE`.

- nsim:

  A single positive integer, the sample size when `approx = "mc"`. Read
  only when `expected = TRUE`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of fifteen numeric vectors named for the multi-index they
carry, from `a_a_a_a` to `p_p_p_p`, each of length
`max(length(y), lengths(theta))`.

## Notation

\\\ell\\ is the log-density of one observation, \\a \> 0\\ the scale,
\\d \> 0\\ and \\p \> 0\\ the two shapes.

## See also

[`distrib_deriv3.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
for the order below and the construction,
[`gengamma_components()`](https://statmodels7.github.io/distributions7/reference/gengamma_components.md)
for the assembly, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
y <- c(0.6, 1.4, 3.1)
th <- list(a = 2, d = 1.5, p = 1.3)
d4 <- distrib_deriv4(d, y, th)
length(d4)
#> [1] 15
names(d4)[1:4]
#> [1] "a_a_a_a" "a_a_a_d" "a_a_a_p" "a_a_d_d"

# A central difference of the third order reproduces the pure-scale
# component.
eps <- 1e-4
up <- distrib_deriv3(d, y, list(a = 2 + eps, d = 1.5, p = 1.3))$a_a_a
dn <- distrib_deriv3(d, y, list(a = 2 - eps, d = 1.5, p = 1.3))$a_a_a
all.equal((up - dn) / (2 * eps), d4$a_a_a_a, tolerance = 1e-5)
#> [1] TRUE
```
