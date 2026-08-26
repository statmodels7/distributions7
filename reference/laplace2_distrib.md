# Laplace Distribution, Location and Rate

Builds the distribution object for the Laplace (double exponential)
family written by its **rate** \\\lambda \> 0\\, with density
\\(\lambda/2)\exp(-\lambda\|y-\mu\|)\\. It is the same law as
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
at \\\lambda = 1/\sigma\\, and every quantity that does not involve the
second parameter is identical between them.

The rate form is the one penalized regression uses: with \\\mu = 0\\ the
negative log-density is \\\lambda\|y\| - \log(\lambda/2)\\, so
\\\lambda\\ is a shrinkage parameter and larger values penalize harder.
`penalties7::lasso_penalty()` is this family with the location held at
zero.

## Usage

``` r
laplace2_distrib(link_mu = identity_link(), link_lambda = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the location \\\mu\\.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the location ranging over the whole line already.

- link_lambda:

  A `link` object from `linkfunctions7` for the rate \\\lambda\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `Laplace2Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"laplace2"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "lambda")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, `link_params` the two links given here, and
`params_smooth` `c(mu = FALSE, lambda = TRUE)`.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \lambda) =
\dfrac{\lambda}{2}\exp\left(-\lambda\|y-\mu\|\right),\$\$ with \\\mu \in
(-\infty, \infty)\\ and \\\lambda \in (0, \infty)\\. The mean and median
are \\\mu\\, the variance is \\2/\lambda^2\\, the skewness is 0 and the
excess kurtosis is 3.

Against
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
the map is \\\lambda = 1/\sigma\\. The location is shared, so its score,
its information and its derivatives in the response are the same
quantities read at the same point; the second parameter's are not, and
this parametrization's are the simpler of the two, \\\lambda\\ being the
natural parameter of the exponential family in \\\|y - \mu\|\\.

## The kink

\\\|y - \mu\|\\ is not differentiable at \\y = \mu\\, so this family is
non-regular in its location exactly as
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
is: the observed second derivative in \\\mu\\ is zero, the information
is \\\lambda^2\\, obtained as the variance of the score, and
`params_smooth` is `c(mu = FALSE, lambda = TRUE)`. The consequences are
set out on
[`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md).

A second consequence is visible at the higher orders. The log-density is
\\\log(\lambda/2) - \lambda\|r\|\\, linear in \\\lambda\\ apart from
\\\log\lambda\\ and piecewise linear in \\\mu\\, so **every** third and
fourth derivative except the pure-\\\lambda\\ ones is exactly zero, and
none of them depends on the data. `expected = TRUE` therefore returns
the same values at those orders.

## Estimation

Both estimates are closed form: \\\hat\mu\\ is the sample median and
\\\hat\lambda\\ the reciprocal of the mean absolute deviation about it.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches them by Fisher scoring, which inverts the information above.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\lambda \> 0\\ the rate, with variance \\2/\lambda^2\\. Here
\\\lambda\\ is a distribution parameter; the same letter names a
smoothing parameter in `penalties7` and above, and the two meet in the
lasso, which is this family with \\\mu\\ held at zero.

## References

Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001). *The Laplace
Distribution and Generalizations*. Birkhauser, Boston.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
*Journal of the Royal Statistical Society, Series B*, **58**(1),
267–288.

## See also

[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
for the scale parametrization;
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md),
which mixes this family with a Gaussian and contains it as the
pure-\\\ell_1\\ case;
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
which holds the location at zero to produce a prior;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[Laplace2Distrib](https://statmodels7.github.io/distributions7/reference/Laplace2Distrib.md)
for the class.

## Examples

``` r
d <- laplace2_distrib()
d
#> Distribution: Laplace2
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu     (location)           | Link: identity   | Domain: (-Inf, Inf)  [non-smooth log-likelihood]
#>   lambda (rate)               | Link: log        | Domain: (0, Inf)

# The same law as the scale parametrization at lambda = 1/sigma.
y <- c(-1.2, 0.3, 2.5)
all.equal(distrib_pdf(d, y, list(mu = 0.4, lambda = 1 / 1.5)),
          distrib_pdf(laplace_distrib(), y, list(mu = 0.4, sigma = 1.5)))
#> [1] TRUE

# The variance is 2/lambda^2: a larger rate is a tighter distribution.
vapply(c(0.5, 1, 2), function(l) variance(d, list(mu = 0, lambda = l)),
       numeric(1))
#> [1] 8.0 2.0 0.5

# Both estimates are closed form.
set.seed(12)
z <- distrib_rng(d, 2000, list(mu = 3, lambda = 0.5))
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = median(z),
                 lambda = 1 / mean(abs(z - median(z)))))
#>              mu    lambda
#> fitted 2.967076 0.5157146
#> closed 2.967076 0.5157146

# Holding the location at zero gives the prior a lasso penalty is written
# from: the negative log-density is lambda |y| up to a constant.
prior <- fixed(d, mu = 0)
prior@params
#> [1] "lambda"
-distrib_pdf(prior, c(-2, 0, 2), list(lambda = 0.5), log = TRUE)
#> [1] 2.386294 1.386294 2.386294
0.5 * abs(c(-2, 0, 2)) - log(0.5 / 2)
#> [1] 2.386294 1.386294 2.386294
```
