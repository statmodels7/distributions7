# Every Way to Write an Integer as an Ordered Sum

The weak compositions of `n` into `k` parts: every vector of `k`
non-negative integers summing to `n`, one per row.

## Usage

``` r
compositions(n, k)
```

## Arguments

- n:

  The total, a non-negative integer.

- k:

  The number of parts, a positive integer.

## Value

An integer matrix with `k` columns.

## Details

Built by recursion on the number of parts, which is what keeps the
result in a fixed order and avoids generating and filtering a full grid.
There are `choose(n + k - 1, k - 1)` of them, so the enumeration is only
practical for a moderate size: at `n = 20` and `k = 5` it is 10626 rows,
and at `k = 10` it is 10015005.

## See also

[`mv_support`](https://statmodels7.github.io/distributions7/reference/mv_support.md)

## Examples

``` r
compositions(3, 2)
#>      [,1] [,2]
#> [1,]    0    3
#> [2,]    1    2
#> [3,]    2    1
#> [4,]    3    0
```
