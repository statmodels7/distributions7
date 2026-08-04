# Skew t Analytical Gradient

First derivatives of the log-density. Three of the four are closed form:
with \\D = A + QB\\ in the notation of
[`skewt_pieces`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md),
\$\$\dfrac{\partial \ell}{\partial \mu} = -\dfrac{D}{\sigma}, \qquad
\dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1 + zD}{\sigma}, \qquad
\dfrac{\partial \ell}{\partial \alpha} = Q z c.\$\$

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

A named list of first derivatives.

## Details

The derivative in \\\nu\\ is **not** closed form. It contains \\\partial
\log T\_{\nu+1}(w)/\partial \nu\\, a derivative of the Student \\t\\
distribution function with respect to its degrees of freedom, which has
no elementary expression — the same obstruction the gamma and beta
distribution functions meet in their shape. That one component is
obtained by a single central difference of the log-density, never by
differencing another difference.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
