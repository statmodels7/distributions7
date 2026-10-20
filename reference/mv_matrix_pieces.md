# The Matrix, Its Inverse and Their Derivatives

Computes the matrix quantities from which the elliptical families
evaluate, for a matrix parametrization that carries either \\\Sigma\\ or
\\\Sigma^{-1}\\: the matrix, its inverse, the log-determinant of
\\\Sigma\\, and the first and second derivative arrays of \\\Sigma\\ and
of \\\Sigma^{-1}\\.

## Usage

``` r
mv_matrix_pieces(s, eta, inverted, derivs = FALSE, derivs2 = FALSE)
```

## Arguments

- s:

  A parameters7 matrix parametrization.

- eta:

  Its free vector.

- inverted:

  `TRUE` when `s` carries \\\Sigma^{-1}\\.

- derivs, derivs2:

  Logical flags that request the first and the second derivative arrays.
  `derivs2 = TRUE` computes the first derivatives as well.

## Value

A list with `sigma`, `sigma_inv`, `logdet` (of \\\Sigma\\) and
`inverted`; with `derivs` or `derivs2` also `a` and `ai`, the first
derivative arrays of \\\Sigma\\ and of \\\Sigma^{-1}\\, and `dlogdet`;
with `derivs2` also `a2`, `ai2` and `d2logdet`.

## Details

Where the parametrization carries \\\Sigma\\, the derivatives of
\\\Sigma\\ are its own arrays and those of \\\Sigma^{-1}\\ come from
[`parameters7::param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.html)
and
[`parameters7::param_inv_d2()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.html),
which are exact for
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
where the matrix is nearly singular. Where it carries \\\Sigma^{-1}\\,
the derivatives of the inverse are the parametrization's own arrays, and
those of \\\Sigma\\ follow from the chain rule for an inverse,
\$\$\partial_k\Sigma = -\Sigma A_k \Sigma, \qquad \partial\_{kl}\Sigma =
\Sigma\left(A_l\Sigma A_k + A_k\Sigma A_l - A\_{kl}\right)\Sigma,\$\$
with \\A_k\\ and \\A\_{kl}\\ the parametrization's arrays.

The derivatives of the log-determinant are returned for \\\Sigma\\ in
either case, as `dlogdet` and `d2logdet`. Where the parametrization
carries \\M = \Sigma^{-1}\\, \\\log\|\Sigma\| = -\log\|M\|\\, so the
parametrization's own derivatives are negated here, once for both
families.

## See also

[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
and
[`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md),
the two callers.
