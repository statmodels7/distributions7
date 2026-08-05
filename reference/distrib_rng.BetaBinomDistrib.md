# Beta-Binomial Random Generation

Draws a probability from the Beta and then a Binomial with it, which is
the hierarchy the family is defined by.

## Arguments

- distrib:

  A `BetaBinomDistrib` object.

- n:

  The number of draws.

- theta:

  A list containing `mu` and `sigma`.

## Value

A numeric vector of length `n`.

## See also

[`betabinom_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom_distrib.md)
