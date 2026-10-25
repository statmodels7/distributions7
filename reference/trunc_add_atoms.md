# The Atoms of a Truncated Continuous Family in an Expectation

Adds to an expectation over the continuous part of a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family the masses of the parent's atoms inside the retained interval,
each multiplied by `f` at the atom.

## Usage

``` r
trunc_add_atoms(distrib, f, theta, cont, ...)
```

## Arguments

- distrib:

  A truncated continuous family.

- f:

  The function whose expectation is taken.

- theta:

  A named list of parameters.

- cont:

  The expectation over the continuous part.

- ...:

  Further arguments passed to `f`.

## Value

A numeric vector with one value per parameter combination.
