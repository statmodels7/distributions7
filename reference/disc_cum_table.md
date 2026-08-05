# Cumulative Probability Table for a Discrete Distribution

Builds the cumulative pmf from the support lower bound, growing the
table geometrically until it covers the probability and the support
point asked for.

## Usage

``` r
disc_cum_table(distrib, theta, need_p = -Inf, need_k = -Inf, kmax = 1e+06)
```

## Arguments

- distrib:

  An object inheriting from class `"discrete_distrib"`.

- theta:

  A named list of parameters.

- need_p:

  Grow the table until it covers at least this probability.

- need_k:

  Grow the table until it reaches at least this support point.

- kmax:

  A ceiling on the table size, to bound a runaway request.

## Value

A list holding the support points and their cumulative probabilities.

## Details

This is the whole of the discrete fallback: the cdf, the quantile
function and the random generator are all lookups into it. No new
algorithm was needed for discrete distributions, because the cdf of a
discrete variable is a step function and inverting it is exact. The cost
was one R-level call per draw, and the lookup is vectorized.

Requires a finite lower bound, which every standard count distribution
has.
