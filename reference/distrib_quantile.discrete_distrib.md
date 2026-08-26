# Default Numerical Quantile Function for Discrete Distributions

The fallback for a discrete family that implements no analytical
quantile function: the smallest support point \\k\\ with \\F(k) \ge p\\,
read off the cumulative table
[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md)
builds. It is **exact**. The distribution function of a lattice variable
is a step function, so inverting it is a lookup and there is nothing to
solve; on a Poisson defined by its mass alone it returns
[`stats::qpois()`](https://rdrr.io/r/stats/Poisson.html)'s answer at
every probability.

Note the consequence a reader coming from the continuous case needs: the
round trip does not close. Asking for 0.025, 0.5 and 0.975 returns
points whose cumulative probabilities are larger, because no support
point sits exactly at those values.

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` that registers no method
  of its own.

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
  under `log.p`.

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
[`distrib_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.discrete_distrib.md),
the function it inverts;
[`distrib_rng.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.discrete_distrib.md),
which reads the same table;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.
