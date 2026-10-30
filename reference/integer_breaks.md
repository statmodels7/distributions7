# Whole-Number Breaks of an Axis

Returns the breaks that [`pretty()`](https://rdrr.io/r/base/pretty.html)
chooses over `limits`, keeping the whole numbers, for the horizontal
axis of a graph over the support of a discrete family.

## Usage

``` r
integer_breaks(limits)
```

## Arguments

- limits:

  The limits of the axis, a numeric vector of length 2.

## Value

A numeric vector of whole numbers.
