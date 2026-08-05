# Beta-Binomial Analytical Observed Hessian

The two-variable chain rule from the shapes, whose second derivatives
are differences of trigammas; the mixed shape component carries only the
\\S = \alpha + \beta\\ part, the two shapes entering the mass function
separately otherwise.

## Arguments

- distrib:

  A `BetaBinomDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `sigma`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list of second-derivative components.

## See also

[`betabinom_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom_distrib.md)
