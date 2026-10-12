# Upper-Tail Table for a Discrete Distribution

Builds \\S(k) = P(Y \> k)\\ over the support from its lower bound by
summing the mass from the top of the table downward, so that each value
is a sum of the terms above \\k\\ and keeps its relative accuracy
however small it is. The table grows by doubling until it reaches the
upper bound of the support or until the mass beyond its last point is
negligible, below \\10^{-17}\\ of the probability that the caller needs.

## Usage

``` r
disc_upper_table(distrib, theta, need_p = NULL, need_k = NULL, kmax = 1e+06)
```

## Arguments

- distrib:

  An object inheriting from class `"discrete_distrib"`.

- theta:

  A named list of parameters, each of length 1.

- need_p:

  The probability that the caller compares with the table: the table
  grows until the neglected mass is below \\10^{-17}\\ of it.

- need_k:

  A support point whose upper tail the caller reads: the table grows
  until it passes this point and the neglected mass is below
  \\10^{-17}\\ of the mass above it.

- kmax:

  A ceiling on the table size.

## Value

A list with the support points `ks` and the upper-tail probabilities
`S`, of the same length.

## Details

The mass beyond the last point \\K\\ is bounded by the geometric series
\\f(K)\\r/(1-r)\\ with \\r = f(K)/f(K-1)\\, which holds where the ratio
of consecutive masses decreases, as it does in the tail of every shipped
count family, and is added to every entry. A ratio that has not fallen
below one keeps the table growing up to `kmax` points, after which the
function warns.

## See also

[`disc_cum_table()`](https://statmodels7.github.io/distributions7/reference/disc_cum_table.md),
the table of the lower tail;
[`distrib_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.discrete_distrib.md)
and
[`distrib_quantile.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.discrete_distrib.md),
which read it.
