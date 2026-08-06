# Derivative Components of the Generalized Gamma

The components of \\\partial^{\alpha+\beta+\gamma}\ell / \partial
a^\alpha \partial d^\beta \partial p^\gamma\\ at any order from one to
four, assembled from the five terms of the log-density.

## Usage

``` r
gengamma_components(y, theta, order)
```

## Arguments

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

- order:

  The derivative order, 1 to 4.

## Value

A named list of component vectors, keyed as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## See also

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md),
[`fdb2`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
