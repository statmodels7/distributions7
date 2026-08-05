# Generalized Gamma Analytical Expected Hessian

Closed form. Every expectation the observed Hessian needs is a moment of
\\u = (Y/a)^{p}\\, which is Gamma with shape \\k = d/p\\ and unit rate:
\\\mathbb{E}\[u\] = k\\, \\\mathbb{E}\[u\log u\] = k\psi(k+1)\\ and
\\\mathbb{E}\[u(\log u)^{2}\] = k\\\psi(k+1)^{2} + \psi'(k+1)\\\\, so
`approx` is ignored.

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

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

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
