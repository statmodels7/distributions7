# Mean of a Distribution Object

Computes the expected value \\E\[Y\]\\ of a distribution object
numerically via
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).
Distribution classes with a closed-form mean may override this method
with an analytical version.

## Arguments

- x:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameters. Vectors are supported.

- ...:

  Additional arguments passed to
  [`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).

## Value

A numeric vector of means.
