# Multinomial Analytical Third and Fourth Derivatives

Closed form. The log-mass is \\\sum_j y_j \log p_j\\ up to a constant,
so each term depends on one coordinate of the simplex and the chain rule
is one univariate partition sum per coordinate, with \\f^{(m)}(p) =
(-1)^{m-1}(m-1)!\\p^{-m}\\.

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- expected:

  Whether to return expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  The approximation used when `expected` is `TRUE`.

- nsim:

  The Monte Carlo sample size.

- ...:

  Unused.

## Value

A named list of third-derivative components.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
