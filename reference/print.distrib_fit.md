# Print Method for Maximum-Likelihood Fits

Shows a fit in four blocks: the family and the sample size with the
log-likelihood and the two information criteria; what the optimizer did;
the estimates on the parameter scale; and the same on the link scale, so
that the interval the fit actually built is visible beside its image.

The optimizer line names the method that produced the estimates, which
is not always the one asked for, and the convergence line names the
**stopping rule** that ended the run. Without that rule `converged` says
nothing: it records that some test was met. A run that did not converge
also prints the score per observation at the point it stopped at, which
is the one number that says whether the point is usable.

For a multivariate fit the coordinates of the covariance structure are
replaced by the quantities the model is written in — the location, the
standard deviations and the correlations that
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
derives — because nobody reads a log-Cholesky coordinate. Any parameter
the structure does not account for, such as the degrees of freedom of a
\\t\\, is printed after them.

## Arguments

- x:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- digits:

  Number of significant digits for every table and for the header
  figures. Defaults to 4. Passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html) for the tables,
  so it is a number of decimal places there.

- ...:

  Unused, accepted for compatibility with
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly. Called for the output it writes.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
for the object;
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md)
to recompute the intervals at another level;
[`plot.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/plot.distrib_fit.md)
to compare the fit with the data;
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for the multivariate table.

## Examples

``` r
set.seed(1)
d <- gaussian1_distrib()
y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
print(fit_distrib(d, y))
#> Maximum-likelihood fit: gaussian1
#> Observations: 500   Log-likelihood: -1264   AIC: 2532   BIC: 2541
#> Method: Fisher scoring   iterations: 2   evaluations: f 3, g 3   time: 2 ms
#> Converged: yes (gradient (max-norm) < 1e-06)
#> 
#> Parameter scale:
#>       Estimate Std. Error   2.5%  97.5%
#> mu      2.0679     0.1356 1.8021 2.3338
#> sigma   3.0327     0.0959 2.8505 3.2267
#> 
#> Link scale:
#>       Estimate Std. Error   2.5%  97.5%
#> mu      2.0679     0.1356 1.8021 2.3338
#> sigma   1.1095     0.0316 1.0475 1.1714

# The link-scale block is where the interval is built. On a Gamma written
# in the mean and the dispersion both parameters carry a log link, so both
# lower limits are positive by construction.
g <- fit_distrib(gamma2_distrib(), rgamma(300, shape = 4, rate = 2))
print(g, digits = 3)
#> Maximum-likelihood fit: gamma2
#> Observations: 300   Log-likelihood: -403   AIC: 810   BIC: 817
#> Method: Fisher scoring   iterations: 3   evaluations: f 4, g 4   time: 3 ms
#> Converged: yes (gradient (max-norm) < 1e-06)
#> 
#> Parameter scale:
#>        Estimate Std. Error  2.5% 97.5%
#> mu        1.948      0.059 1.836 2.067
#> sigma2    1.045      0.103 0.861 1.269
#> 
#> Link scale:
#>        Estimate Std. Error   2.5% 97.5%
#> mu        0.667      0.030  0.607 0.726
#> sigma2    0.044      0.099 -0.149 0.238
```
