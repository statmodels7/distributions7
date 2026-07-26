# Atoms of a Truncated Continuous Distribution

Truncation preserves the parent's atoms that survive it, rescaled by
\\1/Z\\. This matters only when the parent is itself mixed, as
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
of a continuous distribution is.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- theta:

  A named list of the parent's parameters.

## Value

A list with components `y` and `p`.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
