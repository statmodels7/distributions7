# The Pieces a Weibull Evaluates From

Returns the standardized variable \\z = y/\mu\\, its logarithm and the
quantity \\u = z^{\sigma}\\. Every derivative of the Weibull
log-density, in the parameters and in the response, is a polynomial in
\\u\\ and \\u \log z\\, so the three are computed once and shared.

## Usage

``` r
weibull_pieces(y, mu, sigma)
```

## Arguments

- y:

  A numeric vector of observations, each positive. A non-positive value
  gives `NaN` for `lz` and propagates.

- mu:

  The scale parameter, a numeric vector of length 1 or of the length of
  `y`, strictly positive.

- sigma:

  The shape parameter, a numeric vector of length 1 or of the length of
  `y`, strictly positive.

## Value

A named list of three numeric vectors: `z`, the ratio \\y/\mu\\; `lz`,
its logarithm; and `u`, the power \\z^{\sigma}\\. Each has the recycled
length of the inputs.

## Details

\\u\\ is the substitution that makes the family tractable. Under the
model \\u \sim \mathrm{Exp}(1)\\ whatever the parameters are, so an
expectation of any polynomial in \\u\\ and \\\log u\\ is a derivative of
the gamma function at 2, and
[`distrib_expected_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Weibull1Distrib.md)
rests on that. The logarithm is formed before the power so that \\u\\ is
computed as \\\exp(\sigma \log z)\\, which stays representable for a
large shape at a moderate \\z\\.
