# Numerical CDF Derivatives of Any Order

One product stencil of the requested order applied to
[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
which is analytic for every family in the catalog.

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

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- order:

  The derivative order, 3 or 4.

- h_rel:

  Relative finite-difference step.

## Value

A named list of derivative components of \\F\\.

## Details

A repeated parameter contributes the matching higher one-dimensional
factor and distinct parameters each contribute a central two-point
factor, so the whole thing is one stencil rather than a difference of a
difference. The step is \\\varepsilon^{1/(k+2)}\\ scaled by the
parameter, which balances the \\h^{2}\\ truncation against the
\\\varepsilon/h^{k}\\ rounding, and is chosen per observation because
`theta` may vary by observation.

## See also

[`numerical_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md)
