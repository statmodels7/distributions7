# Gaussian Distribution, Mean and Precision

Builds the distribution object for the Gaussian (normal) family
parametrized by its mean \\\mu\\ and its precision \\\tau = 1/\sigma^2
\> 0\\. The returned object carries closed-form derivatives of the
log-density to fourth order, in the parameters and in the response, and
closed-form moments, so every generic of the toolkit answers without a
numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. The defaults are the identity
for the mean, which is already free, and the logarithm for the
precision, which keeps it positive at every predictor.

## Usage

``` r
gaussian3_distrib(link_mu = identity_link(), link_tau = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the mean ranging over the whole line already.

- link_tau:

  A `link` object from `linkfunctions7` for the precision \\\tau\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `Gaussian3Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"gaussian3"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "tau")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \tau) =
\sqrt{\dfrac{\tau}{2\pi}}
\exp\left\\-\dfrac{\tau(y-\mu)^{2}}{2}\right\\,\$\$ with \\\mu \in
(-\infty, \infty)\\ and \\\tau \in (0, \infty)\\. The mean is \\\mu\\,
the variance \\1/\tau\\, and both the skewness and the excess kurtosis
are 0.

This is the same law as
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
in different coordinates, and a separate family for the reason
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
is: the parameter here *is* the precision, and that is what the
estimate, the standard error and the interval describe. The precision is
the parametrization a Bayesian conjugate analysis uses, the gamma being
conjugate for \\\tau\\ at known \\\mu\\.

## Derivatives

Writing \\r = y - \mu\\, the score is \$\$\dfrac{\partial \ell}{\partial
\mu} = \tau r, \qquad \dfrac{\partial \ell}{\partial \tau} =
\dfrac{1}{2\tau} - \dfrac{r^2}{2},\$\$ the observed Hessian
\$\$\ell^{(\mu\mu)} = -\tau, \quad \ell^{(\mu\tau)} = r, \quad
\ell^{(\tau\tau)} = -\dfrac{1}{2\tau^2},\$\$ and its expectation
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\tau, \quad
\mathbb{E}\left\[\ell^{(\mu\tau)}\right\] = 0, \quad
\mathbb{E}\left\[\ell^{(\tau\tau)}\right\] = -\dfrac{1}{2\tau^2}.\$\$
The zero off-diagonal makes the mean and the precision orthogonal, so
their maximum likelihood estimates are asymptotically independent.

It is the flattest of the three parametrizations. The log-density is
quadratic in \\\mu\\ and, in \\\tau\\, is \\\tfrac{1}{2}\log\tau\\ plus
a term linear in \\\tau\\, so every third and fourth derivative is free
of the response: at third order only \\\ell^{(\mu\mu\tau)} = -1\\ and
\\\ell^{(\tau\tau\tau)} = 1/\tau^3\\ survive, and at fourth order only
\\\ell^{(\tau\tau\tau\tau)} = -3/\tau^4\\. Asking those methods for the
expectation returns the same numbers.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale and reaches the
closed-form estimates \\\hat\mu = \bar y\\ and \\\hat\tau = 1 /
\\n^{-1}\sum (y_i - \bar y)^2\\\\, the reciprocal of the maximum
likelihood variance with divisor \\n\\. The example below checks both
against the sample.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\tau = 1/\sigma^2 \> 0\\ the precision. \\\ell^{(ij)}\\ is a second
derivative of \\\ell\\ in parameters \\i\\ and \\j\\. \\\eta\\ is a
parameter on the unconstrained scale of its link, with \\\theta =
g^{-1}(\eta)\\.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. and
Rubin, D. B. (2013). *Bayesian Data Analysis*, 3rd edition, Chapter 2.
Chapman and Hall/CRC, Boca Raton.

## See also

[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
and
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
for the same law in the standard deviation and in the variance;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
for the conjugate prior of \\\tau\\ at known \\\mu\\;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
and
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
for heavier tails;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Gaussian3Distrib](https://statmodels7.github.io/distributions7/reference/Gaussian3Distrib.md)
for the class.

## Examples

``` r
d <- gaussian3_distrib()
d
#> Distribution: Gaussian3
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu  (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   tau (precision)          | Link: log        | Domain: (0, Inf)

# The same law as gaussian1 at sigma = 1/sqrt(tau).
y <- c(-1.2, 0.3, 2.5)
all.equal(distrib_pdf(d, y, list(mu = 1, tau = 0.25)),
          distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#> [1] TRUE

# Moments in closed form: the variance is the reciprocal of the parameter.
th <- list(mu = 1, tau = 0.25)
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#> mean  var skew kurt 
#>    1    4    0    0 

# Fitting recovers the closed-form maximum likelihood estimates.
set.seed(7)
z <- distrib_rng(d, 400, list(mu = 3, tau = 0.25))
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = mean(z), tau = 1 / mean((z - mean(z))^2)))
#>              mu       tau
#> fitted 3.091666 0.2470854
#> closed 3.091666 0.2470854

# The estimate is a precision, so the interval it reports is an interval for
# the precision and stays positive.
confint(fit)
#>          2.5%     97.5%
#> mu  2.8945167 3.2888147
#> tau 0.2151087 0.2838155
```
