# Skew t Analytical Observed Hessian

Second derivatives of the log-density. The block in \\(\mu, \sigma,
\alpha)\\ is closed form; every component involving \\\nu\\ comes from
one finite-difference stencil applied to the log-density, for the reason
given in
[`distrib_gradient.SkewTDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md).

## Arguments

- distrib:

  A `SkewTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma`, `alpha` and `nu`.

- scale:

  Either `"parameter"` or `"link"`.

- ...:

  Unused.

## Value

A named list of second derivatives.

## Details

With \\D = A + QB\\ and \\D' = A' + Q'B^2 + QB'\\, \$\$\dfrac{\partial^2
\ell}{\partial \mu^2} = \dfrac{D'}{\sigma^2}, \qquad \dfrac{\partial^2
\ell}{\partial \mu \\ \partial \sigma} = \dfrac{D + zD'}{\sigma^2},
\qquad \dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1 + 2zD + z^2
D'}{\sigma^2},\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \alpha^2} = Q'
z^2 c^2, \qquad \dfrac{\partial^2 \ell}{\partial \mu \\ \partial \alpha}
= -\dfrac{Q' B z c + Q E}{\sigma}, \qquad \dfrac{\partial^2
\ell}{\partial \sigma \\ \partial \alpha} = -\dfrac{z(Q' B z c + Q
E)}{\sigma}.\$\$ The stencils used for \\\nu\\ are the three-point one
in \\\nu\\ alone and the four-point mixed one otherwise; the mixed
stencil differences two *different* variables, so it is a single stencil
rather than a difference of a difference.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
