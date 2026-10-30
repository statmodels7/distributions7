# The Legend of a Panel Matrix

Builds the legend that
[`mv_pairs_panels()`](https://statmodels7.github.io/distributions7/reference/mv_pairs_panels.md)
places below the panel matrix, as a grob: a graph that draws each
element of the matrix once is built with its legend below the panel, and
the legend is taken from it. The entries are stacked one per row.

## Usage

``` r
mv_pairs_legend(data)
```

## Arguments

- data:

  The data matrix, or `NULL`.

## Value

A gtable holding the legend.
