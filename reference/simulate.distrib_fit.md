# Simulate from a Fitted Distribution

Draws new samples from the fitted distribution, evaluated at the maximum
likelihood estimates. Each replicate has the same length as the data the
model was fitted to, which makes the result directly comparable with the
observations and suitable for a parametric bootstrap or a
posterior-predictive style check.

## Arguments

- object:

  A
  [`distrib_fit`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- nsim:

  Number of replicates to draw. Defaults to 1.

- seed:

  Optional seed. If supplied it is used to initialise the generator, and
  the state of `.Random.seed` in effect before the call is restored
  afterwards, so that simulating does not disturb the calling stream.
  The seed actually used is attached to the result as the `"seed"`
  attribute.

- ...:

  Unused.

## Value

A data frame with `nsim` columns named `sim_1`, ..., `sim_nsim`, each
holding one replicate of `object@n` draws.

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`plot.distrib_fit`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)

## Examples

``` r
set.seed(1)
y <- rnorm(200, 3, 2)
fit <- fit_distrib(gaussian1_distrib(), y)

sims <- simulate(fit, 20, seed = 42)
dim(sims)
#> [1] 200  20

# a parametric bootstrap of any statistic
quantile(vapply(sims, median, numeric(1)), c(0.025, 0.975))
#>     2.5%    97.5% 
#> 2.840640 3.336715 
```
