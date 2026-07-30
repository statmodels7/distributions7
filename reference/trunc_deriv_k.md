# Derivatives of a Truncated Distribution

Builds the order-`k` derivative method for
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md).

## Usage

``` r
trunc_deriv_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method.

## Details

Here \\\ell_T = \ell - \log Z\\, and the ratios are truncated
expectations of the same complete Bell quantity the other wrappers use,
\\d^B Z / Z = \mathbb{E}\_T\[d^B f / f\]\\. Each distinct block costs
one quadrature or summation, which is why they are memoised across the
partition sum.

## See also

[`memo_ratio`](https://statmodels7.github.io/distributions7/reference/memo_ratio.md)
