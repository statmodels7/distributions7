# The Size Behind an NB1 Mean

Returns \\r = \mu/\theta\\, the number of successes that the base R
negative binomial functions take as `size`. Requiring the variance to be
\\\mu(1+\theta)\\ fixes the success probability at \\1/(1+\theta)\\, and
the mean then determines the size. This puts \\\mu\\ inside the gamma
functions of the mass and makes NB1 a different family from the
quadratic-variance one.

The function does not validate its arguments; a non-positive `theta`
gives an infinite or negative size. The methods are reached through the
generics, which signal an error for a non-positive parameter first.

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

A numeric vector of sizes, of the length that the division produces.

## See also

[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for the family, and
[`distrib_pdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin1Distrib.md)
for the mass that it feeds.
