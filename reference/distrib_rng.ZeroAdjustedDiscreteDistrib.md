# Zero-Adjusted Discrete Random Number Generator

Draws zeros with probability \\\pi\\ and otherwise samples from the
zero-truncated parent via inverse transform sampling.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A list with the parent's parameters followed by `za`.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
