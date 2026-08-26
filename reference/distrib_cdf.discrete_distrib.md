# Default Numerical CDF for Discrete Distributions

The fallback for a discrete family that implements no analytical
distribution function: it sums the mass from the finite lower bound of
the support up to \\\lfloor q \rfloor\\. The sum is **exact**, not
approximate, so on a Poisson defined by its mass alone it agrees with
[`stats::ppois()`](https://rdrr.io/r/stats/Poisson.html) to \\10^{-16}\\
absolute, which is the rounding of the addition itself.

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` that registers no method
  of its own.

- q:

  A numeric vector of quantiles; a non-integer is floored, as a
  distribution function on a lattice requires.

- theta:

  A named list of parameters, each of length 1 or `length(q)`.

- lower.tail:

  Logical of length 1. `TRUE`, the default, returns \\P(Y \le q)\\;
  `FALSE` returns the complement.

- log.p:

  Logical of length 1. `TRUE` returns the logarithm, taken after the
  summation.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of cumulative probabilities, the length of the recycled
`q` and `theta`.

## Details

The cumulative table is built once per distinct parameter setting by
[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md)
and every quantile sharing that setting is read from it, so a vector of
quantiles costs one table between them. A support unbounded above is
extended until the remaining mass falls below what the request needs.

## See also

[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md),
which builds the table;
[`distrib_quantile.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.discrete_distrib.md),
which inverts it;
[`distrib_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md),
where the same question needs a quadrature;
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.
