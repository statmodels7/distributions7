# Skew t Probability Density Function

Computes the probability density function, with \\z = (y-\mu)/\sigma\\
and \\w = \alpha z \sqrt{(\nu+1)/(\nu+z^2)}\\: \$\$f(y; \mu, \sigma,
\alpha, \nu) = \dfrac{2}{\sigma}\\ t\_\nu(z)\\T\_{\nu+1}(w)\$\$ with
\\t\_\nu\\ the standard Student \\t\\ density and \\T\_{\nu+1}\\ the
Student \\t\\ distribution function on one more degree of freedom.

## Arguments

- distrib:

  A `SkewTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma`, `alpha` and `nu`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
