# Multivariate Student t Observed Hessian

Closed form, obtained by differentiating the score of
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
once more. Every block picks up a term in \\\partial c/\partial\cdot\\,
because the weight depends on the observation through \\q\\; that
dependence is what distinguishes the family from the gaussian, where
\\c\\ is one and those terms are absent.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)`(distrib@params)`.
