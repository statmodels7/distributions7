# Derivatives of the Expected Information as Exact Sums Over a Finite Support

Sums the moment identities of
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
over the support \\\\0, \dots, N\\\\ against the mass, reading the
family's own observed derivatives.

## Usage

``` r
support_dexpected(distrib, y, theta, order, N, threads = 1L)
```

## Arguments

- distrib:

  A discrete distribution with a finite support.

- y:

  The response, read for its length.

- theta:

  A named list of parameters, each of length one or `length(y)`.

- order:

  `1L` or `2L`.

- N:

  The largest point of the support.

- threads:

  Passed to the derivative generics.

## Value

A named list on the parameter scale, keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
at order 1 and as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
followed by
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
at order 2, so that
[`dexpected_analytic()`](https://statmodels7.github.io/distributions7/reference/dexpected_analytic.md)
reads both orders from one call.

## Details

The support being finite, each sum is the expectation exactly: no
quadrature, no truncation. Every observation's parameters are crossed
with the \\N + 1\\ support points and the derivative generics are called
once on the whole grid.

## See also

[distrib_dexpected_hessian.BetaBinom1Distrib](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.NegBin1Distrib.md)
