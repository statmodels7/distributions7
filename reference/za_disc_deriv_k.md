# Derivatives of a Zero-Adjusted Discrete Distribution

Builds the order-`k` derivative method for the hurdle form of
[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

## Usage

``` r
za_disc_deriv_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method.

## Details

The likelihood separates completely, so every mixed index vanishes at
every order. At \\y \> 0\\ the parameter part is the parent's derivative
minus that of \\\log(1 - f_0)\\, whose ratios are \\d^B(1-f_0)/(1-f_0) =
-(d^B f_0/f_0) f_0/(1-f_0)\\.
