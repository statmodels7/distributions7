# Weibull Log-CDF Derivatives

Closed form at every order from the survival function \\S =
\exp\\-(q/\mu)^{\sigma}\\\\. Writing \\h = \sigma(\log q - \log\mu)\\
the exponent is \\L = -e^{h}\\, so its partial derivatives are
\\-e^{h}\\ times the complete Bell polynomial in the partials of \\h\\,
and those are elementary: \\\partial^{j}h/\partial\mu^{j} =
\sigma(-1)^{j}(j-1)!/\mu^{j}\\, the same without the \\\sigma\\ when one
index names the shape, and zero when two do.

## Arguments

- distrib:

  A `Weibull1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

- ...:

  Unused.

## Value

A named list, one vector per component.

## See also

[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
