# Recycle a Single Point Against Parameters That Vary by Observation

Checks that `x` and the parameters are conformable and returns `x`
repeated to the common length when it has length one and some parameter
is longer.

## Usage

``` r
recycle_point(distrib, x, theta, arg = "y")
```

## Arguments

- distrib:

  A `distrib` object.

- x:

  The response, quantile or probability, a numeric vector.

- theta:

  The aligned parameter list.

- arg:

  The name of `x` in the calling generic, used in the error message.
  Defaults to `"y"`.

## Value

A numeric vector.

## Details

[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md),
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
call it before dispatch, as the derivative generics call
[`check_derivative_args()`](https://statmodels7.github.io/distributions7/reference/check_derivative_args.md),
and the rule is the same: with \\n\\ the largest of the lengths of `x`
and of the parameters, every parameter must have length 1 or \\n\\ and
`x` must have length 1 or \\n\\; otherwise the function signals an
error. A compiled kernel sizes its output by the length of the response,
so a single point evaluated at parameters of length \\n\\ has to be
repeated \\n\\ times to give one value per observation. This matters
wherever a wrapper evaluates its parent at a fixed point, as
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
does for \\f(0)\\. An empty `x` with parameters of length at most one is
returned as it is. A multivariate response is a matrix and is left
alone.

## See also

[`check_derivative_args()`](https://statmodels7.github.io/distributions7/reference/check_derivative_args.md),
the derivative generics' counterpart.
