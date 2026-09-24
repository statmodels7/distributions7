# The von Mises Families' Derivative Tables

Assembles the keyed components for a two-parameter family in which
nothing moves with the first parameter, \\E\_{11} = f(c)\\, \\E\_{22} =
g(c)\\ and \\E\_{12} = 0\\ in the second parameter \\c\\.

## Usage

``` r
vm_dexpected(P, n, order, f1, f2, g1, g2)
```

## Arguments

- P:

  The two parameter names.

- n:

  The length of the result.

- order:

  `1L` or `2L`.

- f1, f2:

  The first and second derivatives of \\f\\.

- g1, g2:

  The first and second derivatives of \\g\\.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).
