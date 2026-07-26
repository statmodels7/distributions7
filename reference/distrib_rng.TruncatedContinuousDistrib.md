# Truncated Random Number Generator (Continuous)

Inverse transform sampling on the parent: \\Y = Q(F(\ell^-) + U Z)\\
with \\U \sim \mathrm{Unif}(0,1)\\. Exact, and unlike rejection sampling
it always terminates in one pass however small \\Z\\ is.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- n:

  Number of observations to generate.

- theta:

  A named list of the parent's parameters.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md)
