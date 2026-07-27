# Zero-Inflated Third Derivatives

At \\y \> 0\\ the likelihood separates; at \\y = 0\\ it is \\\log L_0\\
with \\L_0\\ affine in \\\zeta\\, so the derivatives follow from the
moment-to-cumulant expansion over set partitions.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list with the parent's parameters followed by `zi`.

- expected:

  Logical; if `TRUE`, the expected derivatives.

## Value

A named list of derivative components.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
