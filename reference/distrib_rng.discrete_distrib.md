# Default Random Generation for Discrete Distributions

The fallback for a discrete family that implements no generator of its
own is inverse transform sampling against the cumulative table that
[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md)
builds.

This method uses no rejection scheme. Inverting the cumulative mass
function is exact, because the distribution function of a discrete
variable is a step function, so there is nothing to solve. With scalar
parameters the table is built once and the whole sample is located in it
by one binary search.

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` that registers no method
  of its own.

- n:

  The number of draws, a single non-negative integer.

- theta:

  A named list of parameters; each has length 1 or `n`. A component of
  length `n` gives one draw per parameter setting, taken in order, and
  the table is rebuilt whenever the setting changes from one draw to the
  next.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws, all of which are support points.

## See also

[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md),
which builds the table;
[`distrib_quantile.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.discrete_distrib.md),
the same lookup at given probabilities;
[`distrib_rng.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.continuous_distrib.md),
where an exact inversion is not available;
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.
