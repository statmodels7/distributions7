# Generalized Pareto Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the generalized Pareto
log-density in \\\sigma\\ and \\\xi\\, **in closed form**, by the
construction
[`distrib_deriv3.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md)
describes carried one order further. This is the order at which the
choice between the two routes bites hardest: the Leibniz form's relative
cancellation is \\(\xi z)^{-b}\\, so at \\b = 4\\ it has lost a third of
the answer by \\\xi z = 1.7\times10^{-4}\\ and everything below that.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead: the expected fourth derivatives have no closed form.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector of observations on the support. With
  `expected = TRUE` only its length is read.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`. `sigma` must be strictly
  positive; `xi` may be of either sign, including zero.

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

- threads:

  A single positive integer, how many threads the polynomial kernel of
  the series branch may use. Defaults to `1L`.

## Value

A named list of five numeric vectors, `sigma_sigma_sigma_sigma`,
`sigma_sigma_sigma_xi`, `sigma_sigma_xi_xi`, `sigma_xi_xi_xi` and
`xi_xi_xi_xi`, each of length `max(length(y), lengths(theta))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\sigma \> 0\\ the
scale, \\\xi\\ the shape, \\z = y/\sigma\\ and \\t = 1 + \xi z\\.

## See also

[`distrib_deriv3.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GPDDistrib.md)
for the order below and the construction,
[`gpd_components()`](https://statmodels7.github.io/distributions7/reference/gpd_components.md)
for the two routes and the measured threshold, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
y <- c(0.5, 2, 8)
th <- list(sigma = 1.5, xi = 0.3)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "sigma_sigma_sigma_sigma" "sigma_sigma_sigma_xi"   
#> [3] "sigma_sigma_xi_xi"       "sigma_xi_xi_xi"         
#> [5] "xi_xi_xi_xi"            

# The two routes agree wherever both are accurate: forcing the cut either
# way at xi * z = 0.05, 0.15 and 0.30.
gc <- distributions7:::gpd_components
gap <- function(xz) {
  xi <- xz / (2 / 1.5)
  a <- gc(2, list(sigma = 1.5, xi = xi), 4L, cut = 1e-9)
  b <- gc(2, list(sigma = 1.5, xi = xi), 4L, cut = 10)
  max(abs(unlist(a) - unlist(b)) / pmax(abs(unlist(a)), 1e-10))
}
vapply(c(0.05, 0.15, 0.30), gap, numeric(1))
#> [1] 1.407712e-08 1.953016e-10 2.495273e-11

# And below the cut the Leibniz form loses the answer entirely, which is
# the reason the series branch exists.
leib <- function(x) gc(2, list(sigma = 1.5, xi = x), 4L, cut = 1e-12)$xi_xi_xi_xi
ser <- function(x) gc(2, list(sigma = 1.5, xi = x), 4L, cut = 10)$xi_xi_xi_xi
rbind(xi = c(1e-2, 1e-4, 1e-6),
      leibniz = vapply(c(1e-2, 1e-4, 1e-6), leib, numeric(1)),
      series = vapply(c(1e-2, 1e-4, 1e-6), ser, numeric(1)))
#>              [,1]         [,2]          [,3]
#> xi       0.010000     0.000100  1.000000e-06
#> leibniz -1.156515 35307.473107  8.561630e+14
#> series  -1.156494    -1.263074 -1.264186e+00
```
