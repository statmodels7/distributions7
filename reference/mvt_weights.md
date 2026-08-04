# The Weight a Multivariate t Gives Each Observation

\\(\nu + p)/(\nu + q)\\ at each observation, with the whitened residuals
that go with it.

## Usage

``` r
mvt_weights(y, pc)
```

## Arguments

- y:

  An \\n \times p\\ matrix.

- pc:

  The result of
  [`mvt_pieces`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md).

## Value

A list with `r`, `w`, `q` and `cw`.

## Details

This weight is the whole of the family's robustness. At \\q = 0\\ it is
\\(\nu+p)/\nu\\ and it decays like \\1/q\\, so an observation far from
the centre contributes less to every derivative rather than dragging the
fit towards itself; letting \\\nu \to \infty\\ sends it to one and
recovers the gaussian, where nothing is downweighted.
