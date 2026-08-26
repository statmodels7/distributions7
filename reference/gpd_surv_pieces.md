# The Exponential Survival Pieces of a Generalized Pareto

Returns \\L = \log(1-F)\\ and an evaluator of its partial derivatives in
\\(\sigma, \xi)\\, in the form
[`register_surv_cdf()`](https://statmodels7.github.io/distributions7/reference/register_surv_cdf.md)
wants. Writing \\L = -z\\\Lambda(u)\\ with \\z = q/\sigma\\ and \\u =
\xi z\\ is what keeps every division by the shape out of the expression.

## Usage

``` r
gpd_surv_pieces(distrib, q, theta)
```

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `sigma` (positive) and `xi` (any real
  value), each a numeric vector of length 1 or `n`.

## Value

A list with `Lval` (a numeric vector), `Lderiv` (a function of a block
of parameter names) and `inside` (a logical vector).

## How the partials split

The scale enters \\z\\ as a plain reciprocal, and \\u\\ is bilinear in
the shape and \\z\\, so a block naming the shape twice contributes
nothing to \\u\\. The partials of \\\Lambda(u)\\ follow by Faa di Bruno
over that, and the product with \\z\\ by Leibniz over the scale indices
alone, the shape not entering \\z\\.

## The support

A negative shape bounds the support above, at \\u = -1\\, so the mask is
`q > 0 & u > -1` and is supplied here, the family's fixed bounds being
unable to see it. Past the upper endpoint every derivative of \\F\\ is
exactly zero.

## Notation

\\\sigma \> 0\\ is the scale, \\\xi\\ the shape of either sign, \\z =
q/\sigma\\, \\u = \xi z\\ and \\\Lambda(u) = \log(1+u)/u\\.

## See also

[`gpd_lambda_derivs()`](https://statmodels7.github.io/distributions7/reference/gpd_lambda_derivs.md)
for the univariate function;
[`surv_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/surv_cdf_deriv_k.md)
and
[`register_surv_cdf()`](https://statmodels7.github.io/distributions7/reference/register_surv_cdf.md);
[`distrib_grad_cdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md)
for the family page.
