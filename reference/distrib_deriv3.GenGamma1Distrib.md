# Generalized Gamma Third and Fourth Derivatives

Closed form at both orders, from
[`gengamma_components`](https://statmodels7.github.io/distributions7/reference/gengamma_components.md):
the log-density is elementary apart from \\\lgamma(d/p)\\ and
\\\exp(p\log(y/a))\\, and each of those is a univariate function of a
two-variable map, so the written-out composition covers every component.

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

- expected:

  Logical; if `TRUE`, the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  The approximation used when `expected` is `TRUE`.

- nsim:

  Monte Carlo draws when `approx = "mc"`.

- ...:

  Unused.

## Value

A named list of third-derivative components.

A named list of fourth-derivative components.

## See also

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
