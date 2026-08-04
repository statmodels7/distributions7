# Multivariate Student t Response Hessian

Closed form. With \\w = \Sigma^{-1}(y-\mu)\\, \\q = (y-\mu)^\top w\\ and
\\c = (\nu+p)/(\nu+q)\\, \$\$\dfrac{\partial^2 \ell}{\partial y \\
\partial y^\top} = -c\\\Sigma^{-1} + \dfrac{2c}{\nu+q}\\ w w^\top,\$\$
which depends on the observation through \\c\\ and \\w\\ — unlike the
gaussian's, which is \\-\Sigma^{-1}\\ everywhere — so one matrix is
returned per row.

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

A \\p \times p \times n\\ numeric array.
