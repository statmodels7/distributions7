# Derivatives of the Truncation Constant via the Parent's CDF

Computes \\d^B Z = d^B F(u) - d^B F(\ell^-)\\ from the parent's cdf
derivatives, or `NULL` when that route is not available.

## Usage

``` r
trunc_mass_derivs(distrib, theta, order)
```

## Arguments

- distrib:

  A truncated distribution object.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 or 2.

## Value

A named list of derivative components of \\Z\\, or `NULL`.

## Details

Replaces one quadrature per component with two calls on the parent.
Gated by
[`has_exact_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md),
and the callers fall back to quadrature on `NULL`.
