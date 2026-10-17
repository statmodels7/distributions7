# Derivatives of a Zero-Inflated Distribution

Builds the order-`k` derivative method for
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

## Usage

``` r
zi_deriv_k(order)
```

## Arguments

- order:

  The derivative order, 3, 4 or 5.

## Value

A function suitable for registering as an S7 method.

## Details

At \\y \> 0\\ the likelihood separates, \\\log(1 - \zeta) + \ell(y)\\,
so a mixed index gives zero and the pure ones are immediate. At \\y =
0\\ it is \\\log L_0\\ with \\L_0 = \zeta + (1-\zeta) f_0\\, which is
affine in \\\zeta\\. Any block containing two or more \\\zeta\\'s
therefore contributes nothing, and this keeps the partition sum small.
With \\w_0 = (1-\zeta) f_0 / L_0\\, the ratios are \\w_0 (d^A f_0/f_0)\\
for a block that contains only parent parameters, and \\-f_0 (d^A
f_0/f_0)/L_0\\ for a block carrying one \\\zeta\\.
