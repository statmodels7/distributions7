# von Mises Random Generation

The rejection algorithm of Best and Fisher (1979), which draws from a
wrapped Cauchy envelope and accepts with a probability that does not
involve the Bessel function at all.

## Arguments

- distrib:

  A `VonMisesDistrib` object.

- n:

  The number of draws.

- theta:

  A list containing `mu` and `kappa`.

## Value

A numeric vector of length `n`, in \\\[-\pi, \pi)\\.

## See also

[`vonmises_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises_distrib.md)
