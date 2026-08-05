# The Derivatives of the Map, by Jets

Runs the map on jets rather than on numbers, so that it returns every
partial derivative of the parent's parameters with respect to the new
ones, to fourth order and exactly.

## Usage

``` r
reparam_jets(distrib, theta, n)
```

## Arguments

- distrib:

  A reparametrized distribution.

- theta:

  A named list of the new parameters, already aligned.

- n:

  The number of observations.

## Value

A list with the layout `lay` and `th`, a list over observations of the
parent's parameters as jets.

## Details

The map is written in ordinary R and knows nothing of this: the
arithmetic operators and the mathematical functions dispatch on the jet
class of parameters7, and the derivatives come out of the same
expression that computes the value.

A parameter that varies from observation to observation needs one jet
per observation, so the common case of scalar parameters is detected and
the map is run once.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
