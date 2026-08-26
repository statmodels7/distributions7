# Beta-Binomial Random Generation

Draws `n` independent beta-binomial counts by the two-stage hierarchy
the family is defined by: a success probability from
\\\mathrm{Beta}(\alpha, \beta)\\ with \\\alpha = \mu/\sigma\\ and
\\\beta = (1-\mu)/\sigma\\, then a count from
\\\mathrm{Binomial}(n\_{\mathrm{trials}}, p)\\ at that probability, one
fresh probability per draw. The draws depend on `.Random.seed` in the
usual way and consume two of R's streams per variate.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

- n:

  A single positive integer, the number of draws. Note that the number
  of **trials** is the object's `size` property, not this argument.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of length `n`. A component of length 1 is recycled.
  `mu` must lie in \\(0, 1)\\ and `sigma` be strictly positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` counts in \\\\0, \dots, size\\\\.

## See also

[`distrib_quantile.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BetaBinom1Distrib.md)
for the inverse-transform route,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters back from a sample, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)

# The sample moments recover the mean and the overdispersed variance.
set.seed(1)
z <- distrib_rng(d, 2e5, th)
rbind(sample = c(mean = mean(z), var = var(z)),
      theoretical = c(mean(d, th), variance(d, th)))
#>                mean      var
#> sample      2.99398 8.390626
#> theoretical 3.00000 8.400000

# The counts are bounded by the trial count, which is the object's size and
# not the argument n.
range(z)
#> [1]  0 10
```
