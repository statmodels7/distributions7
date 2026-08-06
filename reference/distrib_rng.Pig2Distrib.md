# Orthogonal Poisson-Inverse Gaussian Random Generation

The exact mixture sampler of
[`pig1`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Pig1Distrib.md),
at the dispersion
[`pig2_sigma`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
implies.

## Arguments

- distrib:

  A `Pig2Distrib` object.

- n:

  The number of draws.

- theta:

  A list containing `mu` and `alpha`.

## Value

A numeric vector of length `n`.

## See also

[`pig2_distrib`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
