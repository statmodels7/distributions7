# Recycle a Single Point Against Parameters That Vary by Observation

Returns `x` repeated to the common length of the parameters when `x` has
length one and some parameter is longer, and `x` unchanged otherwise.

## Usage

``` r
recycle_point(distrib, x, theta)
```

## Arguments

- distrib:

  A `distrib` object.

- x:

  The response, quantile or probability, a numeric vector.

- theta:

  The aligned parameter list.

## Value

A numeric vector.

## Details

[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
call it before dispatch, as
[`check_derivative_args()`](https://statmodels7.github.io/distributions7/reference/check_derivative_args.md)
recycles the response for the derivative generics. A compiled kernel
sizes its output by the length of the response, so without it a single
point evaluated at parameters of length \\n\\ returned ONE value, read
at the first observation's parameters. Measured on eight families –
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md),
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md),
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md),
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md),
[`gengamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma2_distrib.md)
and
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
– and reachable wherever a wrapper evaluates its parent at a fixed
point:
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
reads \\f(0)\\ that way, and with parameters varying by observation its
derivatives were out by a relative 0.5 to 38 on those parents while its
density, evaluated at the full response, was right. A multivariate
response is a matrix and is left alone.

## See also

[`check_derivative_args()`](https://statmodels7.github.io/distributions7/reference/check_derivative_args.md),
the derivative generics' counterpart.
