# Confidence Intervals for a Maximum-Likelihood Fit

Returns Wald intervals for the estimated parameters. They are built
symmetrically on the link scale, \\\hat\eta \pm z\_{1-\alpha/2}\\
\mathrm{se}(\hat\eta)\\, and mapped through \\g^{-1}\\ when the
parameter scale is asked for, so a limit cannot leave the parameter's
domain: a scale has a positive lower limit and a probability stays
inside \\(0,1)\\. The two ends are sorted after mapping, because a link
need not be increasing.

Any level is available from the stored estimate and standard error
without refitting, so a fit computed at 95% answers at 99% for the cost
of one quantile.

## Arguments

- object:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- parm:

  Which parameters to report, by name or by position. Missing, the
  default, reports all of them. A name that is not a parameter of this
  fit, or a position outside the range, signals an error naming the
  argument.

- level:

  Confidence level, a single number in \\(0, 1)\\. Defaults to the level
  the fit was computed at, `object@level`.

- scale:

  Either `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html); anything
  else signals an error. It is the fourth argument, so name it:
  `confint(fit, "sigma", "link")` passes `"link"` as `level` and fails
  inside the quantile.

- ...:

  Unused, accepted for compatibility with
  [`stats::confint()`](https://rdrr.io/r/stats/confint.html).

## Value

A numeric matrix with one row per requested parameter and two columns,
named for the two tail probabilities as percentages (`"2.5%"` and
`"97.5%"` at the default level). Row names are the parameter names. Both
entries are `NA` where the standard error is.

## See also

[`coef.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/coef.distrib_fit.md)
and
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
which take the same `scale`;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
whose `level` sets the default here.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
fit <- fit_distrib(d, distrib_rng(d, 400, list(mu = 1, sigma = 2)))

confint(fit)
#>            2.5%    97.5%
#> mu    0.8863877 1.265967
#> sigma 1.8070076 2.075626
confint(fit, level = 0.99)     # no refit: wider, from the same estimates
#>            0.5%    99.5%
#> mu    0.8267514 1.325603
#> sigma 1.7680868 2.121316
confint(fit, "sigma")
#>           2.5%    97.5%
#> sigma 1.807008 2.075626

# The interval is built on the link scale and mapped back, so the lower
# limit of a scale is exp() of a real number and cannot be negative.
confint(fit, "sigma", scale = "link")
#>            2.5%     97.5%
#> sigma 0.5916722 0.7302626
lo_eta <- confint(fit, "sigma", scale = "link")[1]
all.equal(confint(fit, "sigma")[1], exp(lo_eta), check.attributes = FALSE)
#> [1] TRUE

# A probability near the boundary: 47 successes in 50 trials.
fb <- fit_distrib(bernoulli_distrib(), rep(0:1, c(3, 47)))
confint(fb)                    # inside (0, 1) by construction
#>         2.5%     97.5%
#> mu 0.8298259 0.9805197
```
