# Put CDF Derivatives on the Requested Tail and Scale

Converts derivatives of the distribution function \\F\\ into derivatives
of whichever tail was asked for, on the natural or the logarithmic
scale. Every route to a cdf derivative in this package produces
derivatives of \\F\\ itself, so the `lower.tail` and `log` arguments are
handled once, here, and no method has to implement four cases.

## Usage

``` r
cdf_tail_scale(distrib, Fq, dF1, dF2 = NULL, lower.tail, log)
```

## Arguments

- distrib:

  An object inheriting from `distrib`. Only its `params` are read, to
  name and pair the components.

- Fq:

  The distribution function at the quantile, a numeric vector.

- dF1:

  A named list of first derivatives of \\F\\, one component per
  parameter, in the parameter order.

- dF2:

  A named list of second derivatives of \\F\\, keyed as
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
  or `NULL` (the default) when only the gradient is wanted.

- lower.tail:

  Is the lower tail wanted? A single logical. `TRUE` leaves the signs
  alone and `FALSE` flips every one of them.

- log:

  Are derivatives of the log probability wanted? A single logical.
  `TRUE` divides by the probability, which returns `-Inf` or `NaN` in a
  tail where that probability has underflowed to zero.

## Value

A named list of numeric vectors: one per parameter when `dF2` is `NULL`,
otherwise one per
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
component. The gradient is not returned alongside the Hessian.

## The two conversions

Switching to the upper tail flips the sign, \\S = 1 - F\\ giving
\\\partial^I S = -\partial^I F\\ at every order. Switching to the log
scale divides by the probability, which at second order brings in the
familiar correction \$\$\partial^2 \log P = \frac{\partial^2 P}{P} -
\frac{\partial P}{P}\\\frac{\partial P}{P}.\$\$ The two commute, and
both are applied to whatever was handed in.

## What the caller must supply

`dF1` is always read, at second order as well: the correction above
needs the first derivatives of the same tail. A caller asking for a
Hessian therefore passes both lists, and passes them as derivatives of
\\F\\ on the natural scale, whichever tail the result is wanted on.

## Notation

\\F\\ is the distribution function, \\S = 1 - F\\ the survival function,
\\P\\ whichever of the two was asked for, and \\\partial^I\\ a
derivative with respect to a multi-index of parameters.

## See also

[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
the generics whose methods all end here;
[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
and
[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
for two of the routes that feed it.
