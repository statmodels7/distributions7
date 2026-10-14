# Lower-Tail Logarithms From the Logarithm of the Argument

`weibull_log_lower()` returns \\\log F = \log(1 - e^{-H})\\ from \\\log
H\\, and `gamma_log_lower()` returns the logarithm of the regularized
lower incomplete gamma function \\P(k, w)\\ from \\\log w\\, so that a
cumulative hazard or an argument that underflows still gives a finite
logarithm. Below \\\log H = -30\\ the first is \\\log H - H/2\\; below
\\\log w = -50\\ the second is \\k\log w - \log\Gamma(k + 1) - w\\k/(k +
1)\\, the first terms of the series of \\P\\; above, they are
`log(-expm1(-H))` and
[`stats::pgamma()`](https://rdrr.io/r/stats/GammaDist.html) on the log
scale.

## Usage

``` r
weibull_log_lower(lh)

gamma_log_lower(lw, k)
```

## Arguments

- lh, lw:

  Numeric vectors, the logarithm of the cumulative hazard or of the
  argument (`-Inf` at zero).

- k:

  The shape of the gamma, recycled.

## Value

A numeric vector of log probabilities.

## Examples

``` r
# a hazard of e^-800 and an argument of e^-800, both below the double range
distributions7:::weibull_log_lower(c(-3, -800))
#> [1]   -3.02479 -800.00000
distributions7:::gamma_log_lower(c(-3, -800), 0.7)
#> [1]   -2.024582 -559.904192
```
