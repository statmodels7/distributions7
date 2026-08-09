# Elastic-Net Observed Information

The second derivatives of the log-density. The data term is linear in
the two rates, so every rate block is a second derivative of \\\log Z\\
carried through the bilinear map \\(\lambda,\alpha) \mapsto (a, c)\\,
whose own cross term contributes the \\\lambda\alpha\\ entry.

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

A named list of Hessian components.

## See also

[`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
