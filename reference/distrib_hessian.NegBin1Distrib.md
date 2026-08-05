# NB1 Analytical Observed Hessian

The two-variable chain rule through \\r = \mu/\theta\\, whose second
derivative in \\r\\ is \\\psi'(y+r) - \psi'(r)\\ and whose mixed term is
\\-1/(1+\theta)\\.

## Arguments

- distrib:

  A `NegBin1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `theta`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second-derivative components.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
