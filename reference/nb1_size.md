# The Size Behind an NB1 Mean

Returns \\r = \mu/\theta\\, the number of successes the base R negative
binomial functions take as `size`. Requiring the variance to be
\\\mu(1+\theta)\\ fixes the success probability at \\1/(1+\theta)\\, and
the mean then determines the size. It is this that puts \\\mu\\ inside
the gamma functions of the mass and makes NB1 a different family from
the quadratic-variance one.

Nothing is validated; a non-positive `theta` gives an infinite or
negative size and the caller sees the failure at the base R function.

## Usage

``` r
nb1_size(mu, theta)
```

## Arguments

- mu:

  The mean, a positive numeric vector.

- theta:

  The dispersion, a positive numeric vector. The two are used
  elementwise, so vectors of different lengths recycle in the usual way.

## Value

A numeric vector of sizes, of the length the division produces.

## See also

[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for the family,
[`nb1_prob()`](https://statmodels7.github.io/distributions7/reference/nb1_prob.md)
for the success probability that goes with it, and
[`distrib_pdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin1Distrib.md)
for the mass they feed.
