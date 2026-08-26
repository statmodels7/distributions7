# Beta-Binomial Random Generation in Its Shapes

Draws `n` independent beta-binomial counts by the two-stage hierarchy
the family is defined by: a success probability from
\\\mathrm{Beta}(\alpha, \beta)\\, then a count from
\\\mathrm{Binomial}(n\_{\mathrm{trials}}, p)\\ at that probability, one
fresh probability per draw. The draws depend on `.Random.seed` in the
usual way and consume two of R's streams per variate.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

- n:

  A single positive integer, the number of draws. Note that the number
  of **trials** is the object's `size` property, not this argument.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled.
  Both must be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` counts in \\\\0, \dots, size\\\\.

## See also

[`distrib_rng.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BetaBinom1Distrib.md)
for the same draw in the mean and dispersion,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the shapes back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
th <- list(alpha = 2, beta = 3)

# The sample moments recover the closed forms.
set.seed(2)
z <- distrib_rng(d, 2e5, th)
rbind(sample = c(mean = mean(z), var = var(z)),
      theoretical = c(mean(d, th), variance(d, th)))
#>                mean      var
#> sample      4.00786 6.002098
#> theoretical 4.00000 6.000000

# The counts are bounded by the trial count, which is the object's size and
# not the argument n.
range(z)
#> [1]  0 10
```
