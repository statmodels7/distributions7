# Where Each Pair of Free Values Sits in a Structure's Second Derivatives

A lookup from a pair of free-value positions to the position of the
corresponding component of `param_d2()`.

## Usage

``` r
struct_pair_lookup(s)
```

## Arguments

- s:

  A parameters7 structure.

## Value

A named list of integers, keyed `"k:l"` with \\k \le l\\.

## Details

Built from parameters7's own enumeration rather than by taking a
component key apart, for the reason that package documents: a free value
whose label contains the separator splits into the wrong number of
pieces.
