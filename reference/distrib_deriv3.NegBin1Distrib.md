# NB1 Third and Fourth Derivatives

Closed form at both orders, from
[`negbin1_components`](https://statmodels7.github.io/distributions7/reference/negbin1_components.md).
In the size \\r = \mu/\theta\\ the log-likelihood is \\G(r) +
rB(\theta) + C(\theta)\\, so the only composite piece is
\\G(\mu/\theta)\\ and its mixed derivatives follow a recursion in the
powers of \\r\\ and the order of \\G\\.

## Arguments

- distrib:

  A `NegBin1Distrib` object.

- y:

  A numeric vector of counts.

- theta:

  A list containing `mu` and `theta`.

- expected:

  Logical; if `TRUE`, the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  The approximation used when `expected` is `TRUE`.

- nsim:

  Monte Carlo draws when `approx = "mc"`.

- ...:

  Unused.

## Value

A named list of third-derivative components.

A named list of fourth-derivative components.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
