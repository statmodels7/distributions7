# Laplace Analytical Observed Hessian

Computes the analytical observed Hessian of the Laplace log-density.
Because the log-density is piecewise linear in \\\mu\\, the second
derivative with respect to \\\mu\\ is **zero** almost everywhere:

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = 0, \qquad
\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} =
-\dfrac{\mathrm{sign}(y-\mu)}{\sigma^2}, \qquad \dfrac{\partial^2
\ell}{\partial \sigma^2} = \dfrac{\sigma - 2\|y-\mu\|}{\sigma^3}\$\$

The degeneracy of \\\partial^2 \ell / \partial \mu^2\\ means
Newton-Raphson cannot update \\\mu\\; use
[`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
(Fisher scoring), which supplies the correct information \\1/\sigma^2\\
for \\\mu\\.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A list containing the vectors of second derivatives.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
