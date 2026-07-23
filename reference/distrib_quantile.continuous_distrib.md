# Default Numerical Quantile Function for Continuous Distributions

Fallback method: continuous distributions that do not implement an
analytical quantile function get one by root-finding on
[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
(which may itself be the numerical fallback). Brackets start from an
approximate mode and expand geometrically, with the step scaled by the
density height at the mode.

## Arguments

- distrib:

  An object inheriting from class `"continuous_distrib"`.

- p:

  A numeric vector of probabilities.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\.

- log.p:

  Logical; if `TRUE`, probabilities are given as logs.
