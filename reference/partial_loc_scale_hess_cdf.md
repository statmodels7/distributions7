# CDF Hessian When Only Some Parameters Are Location-Scale

The second-order companion of
[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md),
and the
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
body the Student t, the pseudo-Huber and the skew t share. The three
components in the location and the scale come from
[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)'s
formulas, and every component touching a shape parameter is differenced.

## Usage

``` r
partial_loc_scale_hess_cdf(
  distrib,
  q,
  theta,
  lower.tail = TRUE,
  log = TRUE,
  ...
)
```

## Arguments

- distrib:

  An object inheriting from `distrib` whose first two parameters are a
  location and a scale, and which has at least one more.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, the location first and the scale second.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
in that enumeration's order.

## Details

The gradient is assembled here as well, on the same split, because the
log-scale conversion in
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
reads it. Both differencing calls pass a `which` argument, so no
component that has a closed form costs a cdf evaluation.

For a family with \\p\\ parameters of which two are the location and the
scale, the closed block is 3 components of the \\p(p+1)/2\\: 3 of 6 for
the Student t and the pseudo-Huber, 3 of 10 for the skew t.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\\ell_y = \partial\log
f/\partial y\\.

## See also

[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md)
for the first order;
[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
for the closed formulas;
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for the differenced components.
