# Assemble a Derivative of the Expected Information From the Factors

Applies Leibniz's rule to \\\tilde E\_{ab} = \sum\_{ij} G\_{ij} J\_{ia}
J\_{jb}\\ along every sorted tuple of new coordinates of the given
length: the positions of the tuple are assigned in all \\3^r\\ ways to
the factors \\G\_{ij}\\, \\J\_{ia}\\ and \\J\_{jb}\\.

## Usage

``` r
chain_leibniz(Gtab, ix, Pn, p, r, n)
```

## Arguments

- Gtab:

  Nested list over the parent pairs, each holding the derivatives of
  \\G\_{ij}\\ indexed by the code of the tuple.

- ix:

  The index tables from
  [`reparam_map_index()`](https://statmodels7.github.io/distributions7/reference/reparam_map_index.md).

- Pn:

  The new parameter names.

- p:

  The number of parent parameters.

- r:

  The order, `3L` or `4L`.

- n:

  The number of observations, to which every component is recycled.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)
over `Pn`.
