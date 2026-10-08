# The Skew Normal's CDF Gradient at Zero Skewness

Returns
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the centered skew normal when some skewness is exactly zero. Those
rows take the limit of the gradient, from \\F = \Phi(z) - \gamma_1
(z^2 - 1)\phi(z)/6 + O(\gamma_1^{4/3})\\ with \\z = (q - \mu)/\sigma\\:
the normal's derivatives in \\\mu\\ and \\\sigma\\, and \\-(z^2 -
1)\phi(z)/6\\ in \\\gamma_1\\; the other rows go through the map.

## Usage

``` r
sn2_grad_cdf_zero(distrib, q, theta, zero, lower.tail, log)
```

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A parameter list.

- zero:

  A logical vector, the rows whose skewness is zero.

- lower.tail, log:

  As in
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md).

## Value

A named list of three numeric vectors, `mu`, `sigma` and `gamma1`.
