# The Matrix, Its Inverse and Their Derivatives, From Either Side

The four matrix quantities every elliptical family evaluates from,
computed from a matrix parametrization that carries either \\\Sigma\\ or
\\\Sigma^{-1}\\: the matrix, its inverse, the log-determinant of
\\\Sigma\\, and the derivative arrays of \\\Sigma\\.

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

  Whether the first and second derivative arrays are wanted.

  The log-determinant's derivatives come back transported too, as
  `dlogdet` and `d2logdet`. They are the one place a caller is likely to
  get the sign wrong: \\\log\|\Sigma\| = -\log\|M\|\\, so both are the
  parametrization's own negated, and reading them from here rather than
  from
  [`parameters7::param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.html)
  is what keeps the flip from being written by hand at each of the four
  call sites – which is how the Student t's gradient came to be out by
  14 while its Hessian was right.

## Value

A list with `sigma`, `sigma_inv`, `logdet` (of \\\Sigma\\) and,
according to the flags, `a`, `dlogdet` and then `a2`, `d2logdet`, all of
\\\Sigma\\.

## Details

Every derivative in this package is written in \\\Sigma\\, so a
parametrization carrying the other side is transported once, here, by
the chain rule for an inverse:

\$\$\partial_k\Sigma = -\Sigma A_k \Sigma, \qquad \partial\_{kl}\Sigma =
\Sigma\left(A_l\Sigma A_k + A_k\Sigma A_l - A\_{kl}\right)\Sigma,\$\$

with \\A_k\\ and \\A\_{kl}\\ the parametrization's own arrays. Writing
it once is what keeps
[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
and
[`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md)
the same arithmetic: the two families differ in what they do with these
quantities, not in how the quantities are obtained.

## See also

[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
and
[`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md),
the two callers.
