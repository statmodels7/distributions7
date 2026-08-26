# Derivatives of a Transformed Distribution

Builds the order-`k` derivative method for
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

## Usage

``` r
trans_deriv_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method.

## Details

There is nothing to do. The Jacobian does not depend on \\\theta\\, so
every derivative of the transformed log-likelihood is the parent's
evaluated at \\x = g^{-1}(y)\\. The methods exist only so that the
numerical fallback is never reached.
