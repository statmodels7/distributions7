# Response Derivatives of the Weibull at Zero

Replaces the entries of `out` at `y = 0` by the derivative of order `k`
of \\(\sigma - 1)\log y - b^{-\sigma} y^{\sigma}\\ there, from
[`dy_log()`](https://statmodels7.github.io/distributions7/reference/dy_log.md)
and
[`dy_pow()`](https://statmodels7.github.io/distributions7/reference/dy_log.md).
The compiled response kernels of
[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md)
divide by \\y\\ and give `NaN` at zero, where the derivative exists when
\\\sigma = 1\\ and is infinite otherwise.

## Usage

``` r
weibull3_dy_at_zero(out, y, theta, k)
```

## Arguments

- out:

  The kernel's result, a numeric vector.

- y:

  The response, of the length of `out`.

- theta:

  A list with the mean and the shape, in that order.

- k:

  The derivative order, 1 to 4.

## Value

`out` with its entries at `y = 0` replaced.
