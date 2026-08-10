# Assemble a Dirichlet's Higher Derivatives

Runs the univariate chain rule over the shape vector for every component
of the requested order and adds the two terms that do not pass through
it: the response, which enters \\\alpha_j\\ linearly, and
\\\log\Gamma(\phi)\\.

## Usage

``` r
dirichlet_higher(distrib, y, theta, order)
```

## Arguments

- distrib:

  A `DirichletDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- order:

  The derivative order, 3 or 4.

## Value

A named list of numeric vectors, one per component.
