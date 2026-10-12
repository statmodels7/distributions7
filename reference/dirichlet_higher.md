# Assemble a Dirichlet's Higher Derivatives

Runs the univariate chain rule over the shape vector for every component
of the requested order and adds the terms that do not pass through it:
the response, which enters \\\alpha_j\\ linearly, \\\log\Gamma(\phi)\\,
and on the additive log-ratio chart the derivatives of \\\sum_j
\log\alpha_j\\ (see
[`distrib_deriv3.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.DirichletDistrib.md)).

## Usage

``` r
dirichlet_higher(distrib, y, theta, order, expected = FALSE)
```

## Arguments

- distrib:

  A `DirichletDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 to 4.

- expected:

  Logical of length 1. When `TRUE` the response term is replaced by its
  expectation, through \\\mathrm{E}\log y_j = \psi(\alpha_j) -
  \psi(\phi)\\.

## Value

A named list of numeric vectors of length `nrow(y)`, keyed as
`deriv_names(distrib@params, order)`.
