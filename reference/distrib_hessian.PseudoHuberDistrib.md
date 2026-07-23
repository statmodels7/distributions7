# Pseudo-Huber Analytical Observed Hessian

Computes the analytical observed Hessian of the Pseudo-Huber
log-density. Let \\r = y - \mu\\, \\D = \sqrt{\nu + (r/\sigma)^2}\\,
\\R_1 = K_1'(\sqrt{\nu})/K_1(\sqrt{\nu})\\ and \\R_2 =
K_1''(\sqrt{\nu})/K_1(\sqrt{\nu})\\:

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{\nu}{\sigma^2
D^3}\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \sigma^2} =
\dfrac{\sigma^4 - 3\sigma^2 r^2 D^{-1} + r^4 D^{-3}}{\sigma^6}\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4D^3} +
\dfrac{1}{2\nu^2} + \dfrac{1}{4}\left(\dfrac{R_1}{\nu^{3/2}} +
\dfrac{R_1^2}{\nu} - \dfrac{R_2}{\nu}\right)\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \mu \partial \sigma} = \dfrac{-2\nu\sigma^2 r -
r^3}{\sigma^2(\nu\sigma^2 + r^2)^{3/2}}\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \mu \partial \nu} = -\dfrac{r}{2\sigma^2 D^3}\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu} =
-\dfrac{r^2}{2\sigma^3 D^3}\$\$

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

## Value

A list containing the vectors of second derivatives.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
