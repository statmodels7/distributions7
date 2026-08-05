# Multinomial Analytical Gradient

\$\$\dfrac{\partial\ell}{\partial\eta_k} = \sum_j
\dfrac{y_j}{p_j}A\_{jk}, \qquad A = \dfrac{\partial p}{\partial\eta}\$\$

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list, one component per parameter.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
