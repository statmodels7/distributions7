# Multinomial Analytical Observed Hessian

\$\$\ell^{(\eta_k\eta_l)} = \sum_j\left(\dfrac{y_j}{p_j}B\_{j,kl} -
\dfrac{y_j}{p_j^2}A\_{jk}A\_{jl}\right)\$\$

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

A named list of second-derivative components.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
