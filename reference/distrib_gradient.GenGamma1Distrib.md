# Generalised Gamma Analytical Gradient

With \\w = (y/a)^{p}\\, \\L = \log(y/a)\\ and \\k = d/p\\,
\$\$\dfrac{\partial\ell}{\partial a} = \dfrac{pw - d}{a}, \qquad
\dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
\dfrac{\partial\ell}{\partial p} = \dfrac{1}{p} +
\dfrac{d\\\psi(k)}{p^{2}} - wL\$\$

## Arguments

- distrib:

  A `GenGamma1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `a`, `d` and `p`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list with the `a`, `d` and `p` components.

## See also

[`gengamma1_distrib`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
