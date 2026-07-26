# Atoms of a Zero-Adjusted Continuous Distribution

The single point mass at zero, with probability \\\pi\\. This is what
makes the object a mixed distribution: its density integrates to \\1 -
\pi\\, not to 1.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object.

- theta:

  A list with the parent's parameters followed by `za`.

## Value

A list with `y = 0` and `p = za`.

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
