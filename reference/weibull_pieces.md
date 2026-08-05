# The Pieces a Weibull Evaluates From

Returns the standardized variable \\z = y/\mu\\, its logarithm and the
quantity \\u = z^{\sigma}\\, which every derivative of the Weibull
log-density is a polynomial in.

## Usage

``` r
weibull_pieces(y, mu, sigma)
```

## Arguments

- y:

  A numeric vector of observations.

- mu:

  The scale parameter.

- sigma:

  The shape parameter.

## Value

A list with `z`, `lz` and `u`.

## Details

\\u\\ is the substitution that makes the family tractable: under the
model \\u \sim \mathrm{Exp}(1)\\ whatever the parameters are, so an
expectation of any polynomial in \\u\\ and \\\log u\\ is a derivative of
the gamma function at 2, which is what
[`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
uses.
