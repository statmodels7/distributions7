# Numerical CDF Derivatives of Any Order

Applies one product stencil of the requested order to
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
which is analytic for every family in the catalog. This is the
continuous family's route above second order, and it is a single
stencil: a repeated parameter contributes the matching one-dimensional
higher-order factor and distinct parameters each contribute a central
two-point factor, so a difference of a difference is never taken.

## Usage

``` r
numerical_cdf_deriv_k(
  distrib,
  q,
  theta,
  order,
  h_rel = .Machine$double.eps^(1/(order + 2))
)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale.

- order:

  The derivative order, 3 or 4.

- h_rel:

  The relative step. It is a single number that defaults to
  \\\varepsilon^{1/(\mathrm{order}+2)}\\. A step much smaller than the
  default is worse, because the rounding grows as
  \\h^{-\mathrm{order}}\\.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## The step

The relative step is \\\varepsilon^{1/(k+2)}\\, which is
\\7.4\times10^{-4}\\ at order 3 and \\2.5\times10^{-3}\\ at order 4, and
it is scaled per component by the parameter. It balances the \\h^2\\
truncation against the \\\varepsilon/h^{k}\\ rounding, and it is chosen
per observation because `theta` may vary by observation.

## Accuracy

Against the closed forms of a location-scale family the stencil has a
relative error, normwise over each component, of the order of
\\10^{-5}\\ at order 3 and \\10^{-4}\\ at order 4. The closed forms are
accurate to rounding error, so a family registers one where it can.

## Notation

\\F\\ is the distribution function, \\\theta\\ the parameter on its own
scale, \\h\\ the step and \\\varepsilon\\ the machine epsilon.

## See also

[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for orders 1 and 2;
[`discrete_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md),
the exact route for a discrete family;
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md),
the closed route against which this is measured.
