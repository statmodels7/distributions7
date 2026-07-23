# Zero-Inflated Random Number Generator

Draws from the wrapped distribution, then replaces a
Bernoulli(\\\zeta\\) fraction of the draws with structural zeros.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list with the parent's parameters followed by `zi`.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
