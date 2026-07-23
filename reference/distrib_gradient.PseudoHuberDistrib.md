# Pseudo-Huber Analytical Gradient

Computes the analytical gradient of the Pseudo-Huber log-density. Let
\\r = y - \mu\\ and \\D = \sqrt{\nu + (r/\sigma)^2}\\:

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2 D}\$\$
\$\$\dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma} \left(
\dfrac{r^2}{\sigma^2 D} - 1 \right)\$\$ \$\$\dfrac{\partial
\ell}{\partial \nu} = -\dfrac{1}{2} \left\[ \dfrac{1}{\nu} +
\dfrac{1}{D} + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\\ K_1(\sqrt{\nu})}
\right\]\$\$

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

## Value

A list containing the vectors of first derivatives.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
