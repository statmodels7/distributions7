# The k-th Response Derivative of an Elementary Term

`dy_log` differentiates \\c\log y\\, `dy_log1m` differentiates
\\c\log(1-y)\\, `dy_pow` differentiates \\c\\y^{p}\\, `dy_logaff`
differentiates \\c\log(a + by)\\ and `dy_cos` differentiates
\\c\cos(y-m)\\.

## Usage

``` r
dy_log(c, y, k)

dy_log1m(c, y, k)

dy_pow(c, p, y, k)

dy_logaff(c, a, b, y, k)

dy_cos(c, m, y, k)
```

## Arguments

- c:

  The coefficient.

- y:

  The response.

- k:

  The derivative order.

- p:

  The exponent, for `dy_pow`.

- a, b:

  The affine coefficients, for `dy_logaff`.

- m:

  The location, for `dy_cos`.

## Value

A numeric vector the length of `y`.

## Details

A term linear in \\y\\ is `dy_pow` at \\p = 1\\ and vanishes from the
second order, so it needs no case of its own.
