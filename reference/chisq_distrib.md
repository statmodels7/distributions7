# Chi-Squared Distribution

Builds the distribution object for the chi-squared family on \\(0,
\infty)\\, parametrized by its mean \\\mu \> 0\\, which is the degrees
of freedom. The returned object carries closed-form derivatives of the
log-density to fourth order, in the parameter and in the response, and
closed-form moments, so every generic of the toolkit answers without a
numerical fallback.

The one argument chooses the link that carries the parameter to the
unconstrained scale an optimizer works on, and defaults to the
logarithm.

## Usage

``` r
chisq_distrib(link_mu = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `ChisqDistrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"chisq"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `"mu"`, `n_params` `1`,
`params_bounds` the domain \\(0, \infty)\\, and `link_params` the link
given here.

## The parametrization

The degrees of freedom are treated as a continuous positive parameter,
so the one-parameter family is estimable. The density on \\y \in (0,
\infty)\\ is \$\$f(y; \mu) = \dfrac{y^{\mu/2 - 1}
e^{-y/2}}{2^{\mu/2}\\\Gamma(\mu/2)},\$\$ the mean is \\\mu\\, the
variance \\2\mu\\, the skewness \\2\sqrt{2/\mu}\\ and the excess
kurtosis \\12/\mu\\. The spread is tied to the mean, so the family has
no free scale.

The law is a gamma with shape \\\mu/2\\ and scale 2, and it is **not** a
gamma with a parameter held: this package writes the gamma in \\(\mu,
\sigma^2)\\, and a scale of 2 is the relation \\\sigma^2 = 2\mu\\
between two parameters rather than a value one of them can be fixed at.
At \\\mu = 2\\ it is the exponential with mean 2.

## The response leaves after the first derivative

The family is a one-parameter exponential family in \\\log y\\, so
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{\log y - \log 2 -
\psi(\mu/2)}{2}, \qquad \ell^{(k)} =
-\dfrac{\psi^{(k-2)}(\mu/2)}{2^{k}}, \quad k \ge 2,\$\$ with \\\psi\\
the digamma function. From the second order on nothing involves the
response at all. Two things follow. On the parameter scale the observed
information is exactly the expected information, and the same holds at
third and fourth order, so asking any of those methods for
`expected = TRUE` returns the same numbers. And \\\mathbb{E}\[\log Y\] =
\psi(\mu/2) + \log 2\\ is what gives the score mean zero.

That coincidence does not carry to the scale a fit optimizes on. The
second-order chain rule adds a term
\\h''(\eta)\\\partial\ell/\partial\mu\\ to the link-scale Hessian; the
expected version drops it because the score has mean zero, and a finite
sample does not. Fisher scoring and Newton's method therefore take
different steps here and agree at the optimum, where the summed score
vanishes. The example on
[`distrib_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ChisqDistrib.md)
measures that difference.

The derivatives of the *distribution* function in the parameter have no
elementary form, the derivative of an incomplete gamma in its shape
being hypergeometric, and are taken by finite difference on the analytic
cdf.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. The estimate solves
\\\psi(\hat\mu/2) = \overline{\log y} - \log 2\\ and has no closed form;
the sample mean is the method of moments starting value and lands beside
it.

## Notation

\\\ell\\ is the log-density of one observation and \\\mu \> 0\\ the
mean, which is also the degrees of freedom. \\\psi\\ is the digamma
function and \\\psi^{(m)}\\ its \\m\\th derivative. \\\eta\\ is the
parameter on the unconstrained scale of its link, with \\\mu =
h(\eta)\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions*, Volume 1, 2nd edition, Chapter 18. Wiley, New
York.

## See also

[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
for the two-parameter family this sits inside at \\\sigma^2 = 2\mu\\,
and
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
for its dispersion form;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the case \\\mu = 2\\;
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
for a wider family still;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameter;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[ChisqDistrib](https://statmodels7.github.io/distributions7/reference/ChisqDistrib.md)
for the class.

## Examples

``` r
d <- chisq_distrib()
d
#> Distribution: Chisq
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (mean)               | Link: log        | Domain: (0, Inf)

# The density is stats::dchisq at df = mu.
y <- c(1, 4, 9)
th <- list(mu = 4)
all.equal(distrib_pdf(d, y, th), dchisq(y, df = 4))
#> [1] TRUE

# Moments: the variance is 2 mu, so the spread is not free.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>     mean      var     skew     kurt 
#> 4.000000 8.000000 1.414214 3.000000 
c(2 * sqrt(2 / 4), 12 / 4)
#> [1] 1.414214 3.000000

# It is the gamma at sigma2 = 2 mu, and the exponential at mu = 2.
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(gamma2_distrib(), y, list(mu = 4, sigma2 = 8)))
#> [1] TRUE
all.equal(distrib_pdf(d, y, list(mu = 2)),
          distrib_pdf(exponential_distrib(), y, list(mu = 2)))
#> [1] TRUE

# From the second order on nothing involves the response, so the observed
# and expected information coincide on the parameter scale.
identical(distrib_hessian(d, y, th), distrib_expected_hessian(d, y, th))
#> [1] TRUE

# Fitting recovers the degrees of freedom; the sample mean starts it off.
set.seed(7)
z <- distrib_rng(d, 2000, th)
rbind(fitted = coef(fit_distrib(d, z)), moment = c(mu = mean(z)))
#>              mu
#> fitted 4.038376
#> moment 4.024857
```
