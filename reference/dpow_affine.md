# Derivatives of a Power of an Affine Function

Returns \\\partial^k u^p/\partial\theta^k\\ where \\u\\ is affine in
\\\theta\\ with slope `inner`, which is
\\(\mathrm{inner})^k\\p(p-1)\cdots(p-k+1)\\u^{p-k}\\. The elastic net's
tail arguments are powers of affine functions of its two
hyperparameters, so this supplies their one-variable derivatives to
[`separable_deriv()`](https://statmodels7.github.io/distributions7/reference/separable_deriv.md).

## Usage

``` r
dpow_affine(u, p, k, inner)
```

## Arguments

- u:

  A numeric vector, the affine function evaluated at the parameter.

- p:

  The power, a single number, not necessarily a whole one.

- k:

  The derivative order, a non-negative whole number. Zero returns
  \\u^p\\ itself.

- inner:

  The slope of the affine function, a single number or a numeric vector
  recyclable against `u`.

## Value

A numeric vector the length of `u`.

## See also

[`separable_deriv()`](https://statmodels7.github.io/distributions7/reference/separable_deriv.md),
which combines two of these;
[`distrib_grad_cdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.EnetDistrib.md),
the consumer.

## Examples

``` r
# The second derivative of (2 theta)^3 at theta = 1 is 8 * 3 * 2 * 2 = 24.
distributions7:::dpow_affine(u = 2, p = 3, k = 2, inner = 2)
#> [1] 48
```
