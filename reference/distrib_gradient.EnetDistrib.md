# Elastic-Net Score

With \\z = y-\mu\\, \\\ell^{(\mu)} = a\\\mathrm{sgn}(z) + cz\\ and the
two rate components are the data terms less the derivatives of \\\log
Z\\, which are \\G/\sqrt{c}\\ and \\-(1+xG)/(2c)\\ for \\a\\ and \\c\\;
the chain to \\(\lambda, \alpha)\\ is linear, \\a = \lambda\alpha\\ and
\\c = \lambda(1-\alpha)\\.

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

A named list of score components.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
