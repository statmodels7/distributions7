# The Success Probability Behind a Geometric Mean

Converts a geometric mean into the success probability that the base R
functions [`stats::dgeom()`](https://rdrr.io/r/stats/Geometric.html),
[`stats::pgeom()`](https://rdrr.io/r/stats/Geometric.html),
[`stats::qgeom()`](https://rdrr.io/r/stats/Geometric.html) and
[`stats::rgeom()`](https://rdrr.io/r/stats/Geometric.html) take as their
`prob` argument.

## Usage

``` r
geom_prob(mu)
```

## Arguments

- mu:

  The mean, a positive numeric vector of any length.

## Value

A numeric vector of probabilities in \\(0, 1)\\, of the length of `mu`.

## Details

The mean number of failures before the first success is \\(1-p)/p\\, so
\\p = 1/(1+\mu)\\. Writing it once keeps the four methods that call
stats from each repeating the algebra, and keeps them from drifting
apart if the parametrization ever changes.

The map is decreasing and carries \\(0, \infty)\\ onto \\(0, 1)\\: a
large mean is a small success probability. No validation is performed,
so a non-positive `mu` returns a value outside \\(0, 1)\\ and the caller
sees the failure at the base R function.

## See also

[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
for the family, and
[`distrib_pdf.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GeometricDistrib.md)
for the first of the four callers.
