# Memoize a Ratio Function on Its Block

Caches a ratio function by the canonical key of its block.

## Usage

``` r
memo_ratio(f, params)
```

## Arguments

- f:

  The ratio function to wrap.

- params:

  The parameter names, in declaration order.

## Value

A function with the same signature as `f`, backed by a cache.

## Details

A partition sum over a fourth-order index evaluates the same small
blocks many times. For the truncated wrapper, a block whose cdf
derivatives are not available in closed form costs a quadrature or a
summation, so memoizing across the partition sum evaluates it once per
distinct block and not once per occurrence.

## See also

[`trunc_deriv_k()`](https://statmodels7.github.io/distributions7/reference/trunc_deriv_k.md)
