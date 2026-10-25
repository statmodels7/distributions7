# Derivative Algebra on Index Multisets

Small operations on quantities represented as functions of a sorted
integer vector `S`, a multiset of parameter indices, returning the
derivative of the quantity along `S` (one value per observation).
`zw34_mul()` forms the product by Leibniz's rule over the sub-multisets
of `S`, `zw34_recip()` the reciprocal from the identity \\A\\(1/A) = 1\\
differentiated along `S`, and the remaining functions are the sum, the
difference and the constants. They serve
[`zero_wrapper_d34expected()`](https://statmodels7.github.io/distributions7/reference/zero_wrapper_d34expected.md).

## Usage

``` r
zw34_memo(f)

zw34_subsets(k)

zw34_mul(A, B)

zw34_recip(A)

zw34_add(A, B)

zw34_neg(A)

zw34_const(x)
```

## Arguments

- f:

  A function of `S` to be memoized.

- k:

  A non-negative integer, the length of `S`.

- A, B:

  Quantities, as functions of a sorted integer vector.

- x:

  A numeric vector, the value of a constant quantity.

## Value

`zw34_memo()`, `zw34_mul()`, `zw34_recip()`, `zw34_add()`, `zw34_neg()`
and `zw34_const()` return a function of `S`; `zw34_subsets()` returns a
list of logical vectors of length `k`, one per sub-multiset.
