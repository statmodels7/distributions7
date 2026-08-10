# The Link Scale of a Second-Order Parameter Derivative

Carries a component keyed by parameter pair onto the unconstrained
scale.

## Usage

``` r
theta2_link_scale(distrib, y, theta, second, first)
```

## Arguments

- distrib:

  A distribution object.

- y:

  The response.

- theta:

  The parameters.

- second:

  The parameter-scale second-order components.

- first:

  The parameter-scale first-order components.

## Value

A named list keyed as `second`.

## Details

The chain rule is diagonal, each parameter having its own link, so a
pair \\(i, j)\\ is multiplied by \\h_i' h_j'\\ and a diagonal pair gains
\\h_i''\\ times the first-order component. The response derivatives do
not enter it: a reparametrization of \\\theta\\ leaves them alone.
