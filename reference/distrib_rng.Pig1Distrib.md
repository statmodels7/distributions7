# Poisson-Inverse Gaussian Random Generation

Exact mixture sampling: \\\lambda\\ is drawn from the inverse Gaussian
with mean \\\mu\\ and shape \\\mu/\sigma\\, whose variance
\\\sigma\mu^2\\ is exactly what the mixing construction requires, and
\\Y \mid \lambda\\ from the Poisson.

## Arguments

- distrib:

  A `Pig1Distrib` object.

- n:

  The number of draws.

- theta:

  A list containing `mu` and `sigma`.

## Value

A numeric vector of length `n`.

## See also

[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
