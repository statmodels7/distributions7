# Numerical Third and Fourth Response Derivatives

One central stencil of the requested order applied to
`distrib_pdf(..., log = TRUE)`.

## Usage

``` r
numerical_deriv_y(
  distrib,
  y,
  theta,
  order,
  h_rel = .Machine$double.eps^(1/(order + 2))
)
```

## Arguments

- distrib:

  An object inheriting from class `"continuous_distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 3 or 4.

- h_rel:

  Relative finite-difference step.

## Value

A numeric vector the length of `y`.

## Details

The stencil reaches two steps either side, so the step is clamped to
half of what
[`fd_steps_y`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
allows: a Gamma observation near zero would otherwise be differentiated
at a point outside the support. The relative step is
\\\varepsilon^{1/(k+2)}\\, which balances the \\h^{2}\\ truncation
against the \\\varepsilon/h^{k}\\ rounding.

## See also

[`numerical_hess_y`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md)

## Examples

``` r
numerical_deriv_y(gaussian1_distrib(), c(-1, 0, 1),
                  list(mu = 0, sigma = 1), order = 3)
#> [1] 0 0 0
```
