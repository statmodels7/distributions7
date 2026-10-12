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
  distribution function of a discrete variable requires.

- theta:

  A named list of parameters, each of length 1 or `length(q)`.

- lower.tail:

  Logical of length 1. `TRUE`, the default, returns \\P(Y \le q)\\;
  `FALSE` returns \\P(Y \> q)\\.

- log.p:

  Logical of length 1. `TRUE` returns the logarithm of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of cumulative probabilities, the length of the recycled
`q` and `theta`.

## Details

The cumulative table is built by
[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md)
and rebuilt whenever the parameter setting of a row differs from that of
the row before, so a vector of quantiles at one setting costs one table.
A support unbounded above is extended, by doubling, until the table
reaches the largest requested support point. With `lower.tail = FALSE`
the upper tail \\P(Y \> q)\\ is \\1 - F(q)\\ where \\F(q) \le 1/2\\;
past the median it is read from
[`disc_upper_table()`](https://statmodels7.github.io/distributions7/reference/disc_upper_table.md),
which sums the mass above \\q\\ from the top of its table, so a tail
probability far below machine epsilon keeps its relative accuracy. A
missing `q` gives `NA`.

## See also

[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md),
which builds the table;
[`disc_upper_table()`](https://statmodels7.github.io/distributions7/reference/disc_upper_table.md),
the table of the upper tail;
[`distrib_quantile.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.discrete_distrib.md),
which inverts it;
[`distrib_cdf.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.continuous_distrib.md),
where the same question needs a quadrature;
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.
