# Index Tuples of a Given Width Over a Number of Variables

The same enumeration
[`parameters7::param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.html)
uses, taken over a count rather than over a parameter, so that a
composite index set formed by appending a coordinate can be enumerated
without building an object for it.

## Usage

``` r
tuple_indices_upto(d, order)
```

## Arguments

- d:

  The number of variables.

- order:

  The tuple width, 1 to 4.

## Value

A list of integer vectors.
