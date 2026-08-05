# NB1 Probability Mass Function

The negative binomial mass at size \\\mu/\theta\\ and success
probability \\1/(1+\theta)\\, which is what makes the variance
\\\mu(1+\theta)\\.

## Arguments

- distrib:

  A `NegBin1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu` and `theta`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`negbin1_distrib`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
