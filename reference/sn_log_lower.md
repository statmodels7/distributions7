# The Skew Normal's Lower Tail on the Log Scale

Returns \\\log\\\Phi(z) - 2T(z, \alpha)\\\\ with \\T\\ Owen's function,
the shape clamped to \\10^{50}\\ in absolute value, past which \\T(z,
\alpha)\\ equals its limit \\\Phi(-\|z\|)/2\\ to double precision. A
value above one half is formed through
[`log1p()`](https://rdrr.io/r/base/Log.html) of the other tail.

## Usage

``` r
sn_log_lower(z, alpha)
```

## Arguments

- z:

  A numeric vector of standardized quantiles.

- alpha:

  A numeric vector of shapes.

## Value

A numeric vector of log probabilities.

## Examples

``` r
distributions7:::sn_log_lower(c(-1, 0, 1), 3)
#> [1] -9.7858035 -2.2787086 -0.3816328
```
