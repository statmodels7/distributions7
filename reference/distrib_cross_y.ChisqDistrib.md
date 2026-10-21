# Chi-Squared Mixed Derivative

The mixed derivative is in closed form: the degrees of freedom are the
mean, and \\\ell^{(y)} = (\mu/2 - 1)/y - 1/2\\ gives \\1/(2y)\\.

## Arguments

- distrib:

  A `ChisqDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`.

- scale:

  Handled by the generic after dispatch; this method returns the
  parameter scale.

- ...:

  Unused.

## Value

A named list with the single component `mu`, a numeric vector of length
`length(y)`.

## See also

[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
for the family;
[`distrib_cross_y.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.Gamma1Distrib.md),
which contains it;
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the generic.
