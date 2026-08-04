# Multivariate Student t Score

Closed form. With \\w = \Sigma^{-1}(y-\mu)\\, \\q = (y-\mu)^\top w\\ and
\\c = (\nu+p)/(\nu+q)\\, \$\$\partial\_\mu \ell = c\\w, \qquad
\partial\_{\eta_k}\ell = -\tfrac{1}{2}\partial\_{\eta_k}\log\|\Sigma\| +
\tfrac{c}{2}\\ w^\top A_k w,\$\$ \$\$\partial\_\nu \ell =
\tfrac{1}{2}\left\[ \psi\\\left(\tfrac{\nu+p}{2}\right) -
\psi\\\left(\tfrac{\nu}{2}\right) - \tfrac{p}{\nu} -
\log\\\left(1+\tfrac{q}{\nu}\right) +
\tfrac{(\nu+p)q}{\nu(\nu+q)}\right\].\$\$ The gaussian score is the
limit \\c \to 1\\.

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

A named list with one numeric vector per parameter.
