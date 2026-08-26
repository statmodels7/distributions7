# Numerical CDF Derivatives of Any Order

One product stencil of the requested order applied to
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
which is analytic for every family in the catalog. This is the
continuous family's route above second order, and it is a single
stencil: a repeated parameter contributes the matching one-dimensional
higher-order factor and distinct parameters each contribute a central
two-point factor, so no difference of a difference is ever taken.

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

  The relative step. A single number, defaulting to
  \\\varepsilon^{1/(\mathrm{order}+2)}\\. A step much smaller than the
  default is worse, the rounding growing as \\h^{-\mathrm{order}}\\.

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

## What it delivers

Measured on a Gaussian against the closed forms that family registers,
the relative error is \\1.5\times10^{-5}\\ at order 3 and
\\1.3\times10^{-4}\\ at order 4. That is a long way short of the closed
route and is the reason a family registers one where it can: the loss is
about five digits at order 3 and about four at order 4.

## Notation

\\F\\ is the distribution function, \\\theta\\ the parameter on its own
scale, \\h\\ the step and \\\varepsilon\\ the machine epsilon.

## See also

[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
for orders 1 and 2;
[`discrete_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv_k.md),
the exact route for a discrete family;
[`loc_scale_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv_k.md),
the closed route this is measured against.
