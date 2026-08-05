# Draw the Panel Matrix of a Multivariate Density

The common engine of
[`plot.multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md)
and the multivariate branch of
[`plot.distrib_fit`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md):
one panel per pair of coordinates, the marginal densities on the
diagonal.

## Usage

``` r
mv_pairs_panels(
  d,
  theta,
  which,
  n_grid,
  col_fit,
  data,
  col_data = "#4682B4",
  ...
)
```

## Arguments

- d:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters, already aligned.

- which:

  The coordinates to show, or `NULL` for all.

- n_grid:

  Points per axis.

- col_fit:

  Color of the fitted density.

- data:

  An \\n \times p\\ matrix of observations, or `NULL`.

- col_data:

  Color of the observed summary.

- ...:

  Unused.

## Value

Invisibly `NULL`.

## Details

When `data` is supplied the diagonal also carries a kernel density
estimate of the observed coordinate and the off-diagonal panels carry
the observations, so that the fitted shape and the sample are read
against each other in the same frame. The kernel estimate is the
comparison that does not assume the model: it is what the data say
without the family's help.
