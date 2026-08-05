# Default Numerical CDF for Continuous Distributions

Fallback method: continuous distributions that do not implement an
analytical CDF get one by numerical integration of `distrib_pdf`. An
approximate mode is located first and the integral is taken over the
side of the mode containing \\q\\ (using the complement for the other
side), so that the quadrature nodes concentrate where the probability
mass is. All quantiles are integrated in one batched call to
[`quad_vec`](https://statmodels7.github.io/numericals7/reference/quad_vec.html),
one row per quantile.

## Arguments

- distrib:

  An object inheriting from class `"continuous_distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logs.

## Value

A numeric vector of cumulative probabilities.
