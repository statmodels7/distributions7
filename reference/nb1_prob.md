# The Success Probability Behind an NB1 Dispersion

Returns \\1/(1+\theta)\\, the success probability the base R negative
binomial functions take as `prob`. It is the value that makes the
variance-to-mean ratio \\1+\theta\\ at every mean, so the dispersion
relative to a Poisson does not change with the level of the counts.

Nothing is validated; a `theta` at or below \\-1\\ gives a value outside
\\(0, 1)\\ and the caller sees the failure at the base R function.

## Usage

``` r
nb1_prob(theta)
```

## Arguments

- theta:

  The dispersion, a positive numeric vector.

## Value

A numeric vector in \\(0, 1)\\, of the length of `theta`, falling
towards 0 as the dispersion grows and to 1 as it goes to zero.

## See also

[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for the family,
[`nb1_size()`](https://statmodels7.github.io/distributions7/reference/nb1_size.md)
for the size that goes with it, and
[`distrib_pdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin1Distrib.md)
for the mass they feed.
