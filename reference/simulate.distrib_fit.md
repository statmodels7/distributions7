# Simulate From a Fitted Distribution

Draws new samples from the fitted distribution, evaluated at the maximum
likelihood estimates and ignoring their uncertainty. Each replicate has
the same length as the data the fit was computed from, so a replicate is
directly comparable with the observations: this is the draw a parametric
bootstrap and a posterior-predictive style check both need.

The draws come from
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
so a family with no closed-form generator is sampled by the package's
own fallback and the cost is that fallback's.

## Arguments

- object:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- nsim:

  Number of replicates to draw, a single positive integer. Defaults to
  1.

- seed:

  Optional seed. When supplied it initializes the generator, and the
  `.Random.seed` in effect before the call is restored afterwards, so
  simulating does not disturb the caller's stream. When `NULL`, the
  default, the caller's stream is used and advanced. Either way the seed
  actually in force is attached to the result as its `"seed"` attribute,
  so a run can be reproduced after the fact.

- ...:

  Unused, accepted for compatibility with
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html).

## Value

A data frame of `object@n` rows and `nsim` columns named `sim_1` to
`sim_nsim`, each column one replicate. The `"seed"` attribute carries
the generator state described above.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
for the fit;
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generator;
[`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)
to compare the fit with the data it came from.

## Examples

``` r
set.seed(1)
y <- rnorm(200, 3, 2)
fit <- fit_distrib(gaussian1_distrib(), y)

sims <- simulate(fit, 20, seed = 42)
dim(sims)
#> [1] 200  20

# A parametric bootstrap of any statistic, here the median.
quantile(vapply(sims, median, numeric(1)), c(0.025, 0.975))
#>     2.5%    97.5% 
#> 2.840640 3.336715 

# The seed argument leaves the caller's stream where it found it.
set.seed(7); before <- runif(1)
set.seed(7); invisible(simulate(fit, 2, seed = 42)); after <- runif(1)
all.equal(before, after)
#> [1] TRUE
```
