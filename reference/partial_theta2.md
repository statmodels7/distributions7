# Splice the Closed Location-Scale Pairs Into the Fallback

Takes the differenced components and replaces the three that the
location-scale identity gives in closed form.

## Usage

``` r
partial_theta2(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  `1` or `2`.

## Value

A named list keyed by parameter pair.
