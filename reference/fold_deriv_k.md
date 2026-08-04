# Derivatives of a Folded Distribution

Builds the order-`k` derivative method for
[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md).

## Usage

``` r
fold_deriv_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method.

## Details

The ratios handed to
[`log_deriv`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
are \\d^B L / L = w\\(d^B f(x)/f(x)) + (1-w)\\(d^B f(-x)/f(-x))\\ with
\\w = f(x)/L\\, each term a complete Bell polynomial in the parent's own
log-derivatives evaluated at one of the two preimages. Folding adds no
parameter, so every index is an index of the parent's.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md),
[`fold_ratio`](https://statmodels7.github.io/distributions7/reference/fold_ratio.md)
