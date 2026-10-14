# The Logarithm of minus log(1 - e^p)

Returns \\\log\\-\log(1 - e^{p})\\\\ for \\p \le 0\\ without forming
\\1 - e^{p}\\: as \\\log(-\log(-\mathrm{expm1}(p)))\\ where \\p \> -\log
2\\, and otherwise as \\p + \log\\-\mathrm{log1p}(-u)/u\\\\, \\u =
e^{p}\\, which is \\p\\ itself once \\u\\ underflows. With \\p = \log
S\\ it is the \\\log(-\log F)\\ of the Gumbel quantile; with \\p = \log
F\\ it is the logarithm of the Weibull cumulative hazard at the
quantile.

## Usage

``` r
log_neg_log1mexp(p)
```

## Arguments

- p:

  A numeric vector of log probabilities.

## Value

A numeric vector of the length of `p`.

## Examples

``` r
# at p = -800 the value is -800 to double precision
distributions7:::log_neg_log1mexp(c(-0.1, -50, -800))
#> [1]    0.8553377  -50.0000000 -800.0000000
```
