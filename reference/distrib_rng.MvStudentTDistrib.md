# Multivariate Student t Generator

The scale-mixture representation: \\\mu + L z / \sqrt{g/\nu}\\ with
\\z\\ standard normal, \\g \sim \chi^2\_\nu\\ and \\LL^\top = \Sigma\\.
A t is a gaussian whose precision has been multiplied by a gamma
variate, which is the same fact that makes it robust.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- n:

  The number of observations to draw.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

An \\n \times p\\ numeric matrix.
