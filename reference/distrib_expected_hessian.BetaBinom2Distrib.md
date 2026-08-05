# Beta-Binomial Analytical Expected Hessian in Its Shapes

An **exact finite sum** over the support, the family being discrete on
\\\\0, \dots, n\\\\: the expectation is a weighted sum of at most
\\n+1\\ terms rather than a quadrature or a sample.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object.

- y:

  A numeric vector of counts.

- theta:

  A list with `alpha` and `beta`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is an exact sum.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`betabinom2_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
