# The Exponential Survival Pieces of a Generalized Pareto

Returns \\L = \log(1-F)\\ and an evaluator of its partial derivatives in
\\(\sigma, \xi)\\.

## Usage

``` r
gpd_surv_pieces(distrib, q, theta)
```

## Arguments

- distrib:

  A `GPDDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `sigma` and `xi`.

## Value

A list with `Lval`, `Lderiv` and `inside`.

## Details

\\L = -z\\\Lambda(u)\\ with \\z = q/\sigma\\ and \\u = \xi z\\. The
scale enters \\z\\ as a plain reciprocal and \\u\\ is bilinear in the
shape and \\z\\, so a block naming the shape twice contributes nothing
to \\u\\; the partials of \\\Lambda(u)\\ follow by Faa di Bruno over
that, and the product with \\z\\ by Leibniz over the scale indices
alone, the shape not entering \\z\\.

## See also

[`gpd_lambda_derivs`](https://statmodels7.github.io/distributions7/reference/gpd_lambda_derivs.md),
[`register_surv_cdf`](https://statmodels7.github.io/distributions7/reference/register_surv_cdf.md)
