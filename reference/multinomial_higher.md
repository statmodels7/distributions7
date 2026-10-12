# Assemble a Multinomial's Higher Derivatives

Computes every component of the requested order. On the additive
log-ratio chart a component of order one is \\y_k - N p_k\\ and a
component of higher order is \\-N\\ times a cumulant from
[`mn_cumulant()`](https://statmodels7.github.io/distributions7/reference/mn_cumulant.md),
with \\N\\ the total of the row. On another chart the method runs the
univariate chain rule over the probability vector for every component,
weighted by the counts of each category.

## Usage

``` r
multinomial_higher(distrib, y, theta, order, expected = FALSE)
```

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 to 4.

- expected:

  Logical of length 1. When `TRUE` each count is replaced by its
  expectation \\n p_j\\, which gives the expected derivatives exactly.

## Value

A named list of numeric vectors of length `nrow(y)`, keyed as
`deriv_names(distrib@params, order)`.
