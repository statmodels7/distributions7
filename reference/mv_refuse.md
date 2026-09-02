# Reject a Quantity That Has No Multivariate Counterpart

Raises the error a multivariate distribution gives for a one-dimensional
quantity, in one wording: the name of the quantity, the family it was
asked of, and a sentence saying what is missing. Every refusal on this
class goes through it, so the eight of them read alike and none of them
silently returns a number of the wrong shape.

## Usage

``` r
mv_refuse(distrib, what, why)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object, whose `distrib_name` is quoted in the message.

- what:

  The name of the quantity, a single string WITHOUT parentheses:
  `"distrib_cdf"`, `"skewness"`. The `()` is appended here, so passing
  them gives `distrib_cdf()()`.

- why:

  One sentence saying what is missing, a single string with no leading
  capital. It is placed after a colon and carries its own final period.

## Value

Never returns: it always signals an error, with `call. = FALSE` so the
message is not prefixed by the internal call.

## See also

[`distrib_cdf.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.multivariate_distrib.md),
[`distrib_quantile.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.multivariate_distrib.md),
[`skewness.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.multivariate_distrib.md)
and
[`kurtosis.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.multivariate_distrib.md)
for the refusals it produces.

## Examples

``` r
d <- mvgaussian1_distrib(2)

# The wording every refusal on this class shares. The name carries no
# parentheses; they are appended here.
try(distributions7:::mv_refuse(d, "some_quantity",
                               "there is no such thing in p dimensions."))
#> Error : some_quantity() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': there is no such thing in p dimensions.

# Which is what a caller sees from the generics themselves.
try(distrib_cdf(d, rbind(c(0, 0)),
                list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
                     sigma_log_L2 = 0, sigma_L2.1 = 0)))
#> Error : distrib_cdf() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': the distribution function is an integral over an orthant, with no closed form and no one-dimensional fallback.
```
