# Split a Multi-Index Into Parent and Wrapper Parts

Separates a multi-index into the parent's parameters and the number of
times the wrapper's own parameter appears.

## Usage

``` r
split_index(idx, mix_name)
```

## Arguments

- idx:

  A character vector of parameter names, with repetition.

- mix_name:

  The wrapper parameter's name.

## Value

A list with `theta`, the parent part, and `n_mix`, a count.

## Details

The wrapper's parameter is always the last one it declares, so the split
needs only its name.
