# Derivatives of a Zero-Adjusted Continuous Distribution

Builds the order-`k` derivative method for the mixed form of
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

## Usage

``` r
za_cont_deriv_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method.

## Details

A continuous parent puts no mass at zero, so there is no truncation
constant to correct for: the parameter part is simply the parent's
derivative, switched off at the atom.
