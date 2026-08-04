# Mean of a Multivariate Student t

The location vector, which is the mean when \\\nu \> 1\\ and undefined
otherwise; `NaN` is returned there rather than the location, since the
location exists as a parameter while the moment does not.

## Arguments

- x:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

A numeric vector of length \\p\\.
