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

A partition of a fourth-order index asks for the same small blocks many
times over. For the truncated wrapper each distinct block costs a
quadrature, so memoizing across the partition sum is the difference
between one integration per block and one per occurrence.

## See also

[`trunc_deriv_k()`](https://statmodels7.github.io/distributions7/reference/trunc_deriv_k.md)
