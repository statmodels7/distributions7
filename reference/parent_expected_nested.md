# The Parent's Derivatives of the Expected Information as Nested Lists

Reads the keyed derivatives of the parent's expected information into
lists indexed by the parent's parameter positions, so that the chain
formulas address \\E\_{ij,k\ldots}\\ by integers.

## Usage

``` r
parent_expected_nested(Pp, src, r)
```

## Arguments

- Pp:

  The parent's parameter names.

- src:

  A named list keyed as
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
  [`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md),
  [`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
  or
  [`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

- r:

  The number of derivative indices, 0 to 4.

## Value

A list nested `2 + r` deep over `seq_along(Pp)`.
