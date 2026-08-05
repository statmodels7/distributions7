# Residuals and Whitened Residuals

The centered response and its image under the inverse covariance, which
are what every derivative of a multivariate gaussian is written in.

## Usage

``` r
mvg_residuals(y, pc)
```

## Arguments

- y:

  An \\n \times p\\ matrix.

- pc:

  The result of
  [`mvg_pieces`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md).

## Value

A list with `r`, the residuals, and `w`, the rows of \\R \Sigma^{-1}\\.
