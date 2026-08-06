# Orthogonal Poisson-Inverse Gaussian Probability Mass Function

The same law as
[`pig1's`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig1Distrib.md)
at \\\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2\\: the parameter
\\\alpha\\ of this parametrization is exactly the argument \\\sqrt{1 +
2\sigma\mu}/\sigma\\ the Bessel function is evaluated at.

## Arguments

- distrib:

  A `Pig2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `alpha`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`pig2_distrib`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
