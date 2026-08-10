# Assemble a Multinomial's Higher Derivatives

Runs the univariate chain rule over the probability vector for every
component of the requested order, weighted by the summed counts.

## Usage

``` r
multinomial_higher(distrib, y, theta, order)
```

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- order:

  The derivative order, 3 or 4.

## Value

A named list of numeric vectors, one per component.
