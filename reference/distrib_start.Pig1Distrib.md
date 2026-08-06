# Poisson-Inverse Gaussian Starting Values

Method of moments: the sample mean for \\\mu\\, and \\(s^2 - \bar
y)/\bar y^2\\ for \\\sigma\\, floored just above zero when the sample is
underdispersed.

## Arguments

- distrib:

  A `Pig1Distrib` object.

- y:

  A numeric vector of observations.

- n_start:

  Ignored; one moment start is returned.

- ...:

  Unused.

## Value

A list with one named parameter list.

## See also

[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
