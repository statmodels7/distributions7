# Multinomial Analytical Expected Hessian

Closed form. \\\mathbb{E}\[y_j\] = n p_j\\, so the first term becomes
\\n\sum_j B\_{j,kl}\\, which vanishes because the probabilities sum to a
constant, and \$\$\mathbb{E}\[\ell^{(\eta_k\eta_l)}\] = -n\sum_j
\dfrac{A\_{jk}A\_{jl}}{p_j}\$\$

## Arguments

- distrib:

  A `MultinomialDistrib` object.

- y:

  A matrix with one row per observation.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Ignored; the expectation is closed form.

- nsim:

  Ignored.

- ...:

  Unused.

## Value

A named list of expected second-derivative components.

## See also

[`multinomial_distrib`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
