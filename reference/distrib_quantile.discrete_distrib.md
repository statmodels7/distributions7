# Default Numerical Quantile Function for Discrete Distributions

The fallback for a discrete family that implements no analytical
quantile function: the smallest support point \\k\\ with \\F(k) \ge p\\,
read off the cumulative table that
[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md)
builds. It is exact. The distribution function of a discrete variable is
a step function, so inverting it is a lookup and there is nothing to
solve. The probability is compared with the running sum after a relative
reduction of \\64\varepsilon\\, as
[`stats::qpois()`](https://rdrr.io/r/stats/Poisson.html) does, so that a
probability equal to \\F(k)\\ returns \\k\\; on a Poisson defined by its
mass alone the result agrees with
[`stats::qpois()`](https://rdrr.io/r/stats/Poisson.html). With
`lower.tail = FALSE` an upper-tail probability below one half is
inverted against the table of \\P(Y \> k)\\ that
[`disc_upper_table()`](https://statmodels7.github.io/distributions7/reference/disc_upper_table.md)
sums from the top, so that it is not turned into \\1 - p\\, which would
lose its digits; the quantile of an upper-tail probability far below
machine epsilon is therefore finite and accurate.

Because the distribution function is a step function, the round trip
does not close. Asking for 0.025, 0.5 and 0.975 returns points whose
cumulative probabilities are larger, because no support point sits
exactly at those values.

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` that registers no method
  of its own.

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  under `log.p`. A value outside \\\[0, 1\]\\ gives `NaN` and a missing
  value `NA`.

- theta:

  A named list of parameters, each of length 1 or `length(p)`.

- lower.tail:

  Logical of length 1. `TRUE`, the default, treats `p` as \\P(Y \le
  q)\\; `FALSE` treats it as the upper tail.

- log.p:

  Logical of length 1. `TRUE` treats `p` as a logarithm and
  exponentiates it first.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of support points, the length of the recycled `p` and
`theta`.

## See also

[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md),
which builds the table;
[`disc_upper_table()`](https://statmodels7.github.io/distributions7/reference/disc_upper_table.md),
the table of the upper tail;
[`distrib_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.discrete_distrib.md),
the function it inverts;
[`distrib_rng.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.discrete_distrib.md),
which reads the same table;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.
