# Append One Block of Derived Quantities to Another

Concatenates two lists in the shape of
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md),
returning the first unchanged when the second is `NULL`.

## Usage

``` r
mv_append_block(out, extra)
```

## Arguments

- out:

  A list as described in
  [`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md).

- extra:

  A list of the same shape, or `NULL`.

## Value

A list as described in
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md).

## See also

[`mv_param_block`](https://statmodels7.github.io/distributions7/reference/mv_derived.multivariate_distrib.md)
