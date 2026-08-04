# Project a Matrix onto What a Structure Can Represent

Returns the free vector of the structure whose matrix is closest to the
one supplied, or the structure's own inverse map when it has one.

## Usage

``` r
struct_free_or_fit(s, m)
```

## Arguments

- s:

  A parameters7 structure.

- m:

  The matrix to represent.

## Value

A numeric vector of length `s@n_free`.

## Details

[`param_free`](https://statmodels7.github.io/parameters7/reference/param_free.html)
is exact or refused: a structure that cannot represent the matrix says
so rather than returning something plausible. That is the right contract
for reporting an estimate and the wrong one for choosing where to begin,
so a refusal here falls back to a short numerical search over the free
values, which is allowed to be approximate because a starting value is.
