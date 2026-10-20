# Assemble a Multivariate Student t's Higher Derivatives

Splits a component into its mean-and-matrix part and its \\\nu\\ part.
The log-determinant contributes only to a component whose indices all
name the matrix, and the orientation of the chart sets its sign (see
[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md)).
A component whose indices all name \\\nu\\ is computed by
[`mvt_cp_d()`](https://statmodels7.github.io/distributions7/reference/mvt_cp_d.md)
and
[`mvt_h_d()`](https://statmodels7.github.io/distributions7/reference/mvt_h_d.md),
which avoid the cancellation of the direct forms at large \\\nu\\. Every
other component is a univariate chain rule in \\q\\ over the set
partitions of its mean and matrix indices, whose outer derivatives
\\(-1)^{m-1}(m-1)!/(\nu+q)^{m}\\ are differentiated in \\\nu\\ by the
Leibniz rule against the linear prefactor \\(\nu+p)/2\\. The two terms
of the Leibniz rule are combined into one fraction, since written
separately they agree to leading order as \\\nu\\ grows.

## Usage

``` r
mvt_higher(distrib, y, theta, order)
```

## Arguments

- distrib:

  A `MvStudentTDistrib` object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 to 4.

## Value

A named list of numeric vectors of length \\n\\, keyed as
`deriv_names(distrib@params, order)`.
