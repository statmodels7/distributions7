# Is This Distribution Already Truncated?

`TRUE` for either of the two truncated classes.

## Usage

``` r
is_truncated(distrib)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

## Value

A single logical.

## Details

Used by
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
to collapse nested truncation to the intersection of the intervals.
Nesting would be correct, but it pays the quadrature cost twice for a
law that a single truncation already describes.
