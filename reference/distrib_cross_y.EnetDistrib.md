# Elastic-Net Mixed Response-Parameter Derivatives

\\\partial^2\ell/\partial y\\\partial\theta_i\\, obtained by
differentiating \\\ell^{(y)} = -a\\\mathrm{sgn}(z) - cz\\ in each
parameter.

## Arguments

- distrib:

  An `EnetDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- scale:

  Either `"parameter"` or `"link"`.

- ...:

  Ignored.

## Value

A named list, one component per parameter.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
