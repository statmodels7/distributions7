# Laplace Distribution, Location and Scale

Builds the distribution object for the Laplace (double exponential)
family with location \\\mu\\ and scale \\\sigma \> 0\\. The returned
object carries closed-form derivatives of the log-density to fourth
order and closed-form moments.

The family is symmetric about \\\mu\\ with variance \\2\sigma^2\\ and
excess kurtosis 3, so it is heavier tailed than a Gaussian. Its maximum
likelihood estimates are the sample median and the mean absolute
deviation about it, both available in closed form.

This family is **not regular** in \\\mu\\: the density has a corner at
\\y = \mu\\ and the observed second derivative there is zero almost
everywhere, while the information is \\1/\sigma^2\\. The object records
this in `params_smooth`, and the details below say what follows from it.

## Usage

``` r
laplace_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the location \\\mu\\.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the location ranging over the whole line already.

- link_sigma:

  A `link` object from `linkfunctions7` for the scale \\\sigma\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `LaplaceDistrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"laplace"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, `link_params` the two links given here, and
`params_smooth` `c(mu = FALSE, sigma = TRUE)`.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \sigma) =
\dfrac{1}{2\sigma}\exp\left(-\dfrac{\|y-\mu\|}{\sigma}\right),\$\$ with
\\\mu \in (-\infty, \infty)\\ and \\\sigma \in (0, \infty)\\. The
distribution function is one exponential on each side of \\\mu\\ and the
quantile function inverts it in closed form.

The mean and the median are \\\mu\\, the variance is \\2\sigma^2\\, the
skewness is 0 and the excess kurtosis is 3.
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
is the same law written by its rate \\\lambda = 1/\sigma\\, which is the
form a lasso penalty uses.

## The kink, and what it costs

\\\|y - \mu\|\\ is not differentiable at \\y = \mu\\, so the log-density
is piecewise linear in \\\mu\\ with a corner. Three consequences, all
visible in the methods:

- the score in \\\mu\\ is \\\mathrm{sign}(r)/\sigma\\, carrying only the
  sign of the residual;

- the observed second derivative in \\\mu\\ is 0 wherever it exists, so
  [`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md)
  returns a vector of zeros;

- the second Bartlett identity fails, and the information is **defined**
  as the variance of the score, \\1/\sigma^2\\. That is what
  [`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md)
  returns.

`params_smooth` is `c(mu = FALSE, sigma = TRUE)`, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
reads it to skip the finite-difference comparison in \\\mu\\, where a
central difference straddling the corner returns a number that is not a
derivative of anything.

## Estimation

Both estimates are closed form: \\\hat\mu\\ is the sample median and
\\\hat\sigma\\ the mean absolute deviation about it, \\n^{-1}\sum\|y_i -
\hat\mu\|\\.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches them by Fisher scoring, which inverts the information above; a
Newton step could not, the observed Hessian being singular in \\\mu\\.
The example below checks the fitted values against both closed forms.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\sigma \> 0\\ the scale, with variance \\2\sigma^2\\. \\r = y - \mu\\
is the residual. The **observed information** is
\\-\partial^2\ell/\partial\theta\\\partial\theta^\top\\ at the data; the
**expected information** is the variance of the score, which for a
regular family equals the expectation of the observed one and here does
not.

## References

Kotz, S., Kozubowski, T. J. and Podgorski, K. (2001). *The Laplace
Distribution and Generalizations*. Birkhauser, Boston.

## See also

[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
for the rate parametrization;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the light-tailed comparison;
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
and
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
for the other robust location families, the second of which smooths this
one's corner away;
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md),
which contains this family as its pure-\\\ell_1\\ case;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
and
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
to estimate the parameters;
[LaplaceDistrib](https://statmodels7.github.io/distributions7/reference/LaplaceDistrib.md)
for the class.

## Examples

``` r
d <- laplace_distrib()
d
#> Distribution: Laplace
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (location)           | Link: identity   | Domain: (-Inf, Inf)  [non-smooth log-likelihood]
#>   sigma (scale)              | Link: log        | Domain: (0, Inf)

# The density, written out.
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
all.equal(distrib_pdf(d, y, th), exp(-abs(y - 0.4) / 1.5) / (2 * 1.5))
#> [1] TRUE

# sigma is a scale: the variance is 2 sigma^2 and the excess kurtosis is 3.
c(variance = variance(d, th), kurtosis = kurtosis(d, th))
#> variance kurtosis 
#>      4.5      3.0 

# Both estimates are closed form: the median and the mean absolute
# deviation about it.
set.seed(12)
z <- distrib_rng(d, 2000, list(mu = 3, sigma = 2))
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = median(z), sigma = mean(abs(z - median(z)))))
#>              mu    sigma
#> fitted 2.967076 1.939057
#> closed 2.967076 1.939057

# The family is not regular in mu: the observed curvature there is 0 and
# the information is 1/sigma^2.
c(observed = mean(distrib_hessian(d, z, th)$mu_mu),
  expected = distrib_expected_hessian(d, 0, th)$mu_mu)
#>   observed   expected 
#>  0.0000000 -0.4444444 
```
