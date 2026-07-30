# Summed Hessian on the Link Scale, as a Matrix

Assembles the package's named list of Hessian components into the
symmetric \\p \times p\\ matrix an optimizer wants, summed over
observations.

## Usage

``` r
fit_hess_matrix(distrib, y, theta, expected)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters on the natural scale.

- expected:

  Logical; whether to use the expected Hessian.

## Value

A symmetric numeric matrix with dimnames taken from the parameters.

## Details

Derivative components are stored one per unique index pair, since the
Hessian is symmetric; this fills both triangles. Passing
`expected = TRUE` gives the expected Hessian, which is what makes Fisher
scoring possible – and what allows a fit on a non-regular family such as
the Laplace, where the observed Hessian is degenerate but the
information is not.

## See also

[`fit_score`](https://statmodels7.github.io/distributions7/reference/fit_score.md),
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
