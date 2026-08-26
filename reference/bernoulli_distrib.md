# Bernoulli Distribution

Builds the distribution object for the Bernoulli family parametrized by
its success probability \\\mu \in (0, 1)\\. The returned object carries
closed-form derivatives of the log-mass to fourth order, observed and
expected, and closed-form moments.

This is the response distribution of logistic regression, and the
default link is the logit that gives that model its name. The logit is
the canonical link here, so the observed and the expected information
coincide on its scale and iteratively reweighted least squares is Fisher
scoring and Newton's method at once.

## Usage

``` r
bernoulli_distrib(link_mu = logit_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the probability \\\mu\\.
  Defaults to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  the canonical link, which maps \\(0, 1)\\ onto the whole line.
  [`linkfunctions7::probit_link()`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.html)
  and
  [`linkfunctions7::cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.html)
  are the usual alternatives; under either the observed and expected
  information differ.

## Value

An S7 object of class `BernoulliDistrib`, inheriting from
`discrete_distrib`, with `distrib_name` `"bernoulli"`, `dimension`
`"univariate"`, `bounds` `c(0, 1)`, `params` `"mu"`,
`params_interpretation` `c(mu = "probability")`, `n_params` `1`,
`params_bounds` the list of \\(0, 1)\\, and `link_params` the one link
given here.

## The parametrization

The mass on \\y \in \\0, 1\\\\ is \$\$P(Y = y; \mu) =
\mu^{y}(1-\mu)^{1-y},\$\$ with \\\mu \in (0, 1)\\. The mean is \\\mu\\,
the variance \\\mu(1-\mu)\\, the skewness \\(1-2\mu)/\sqrt{\mu(1-\mu)}\\
and the excess kurtosis \\(1 - 6\mu(1-\mu))/(\mu(1-\mu))\\. The family
has no dispersion parameter: fixing the mean fixes the variance, which
is why overdispersion in binary data has to come from elsewhere, most
often from a random effect or from grouping.

## The canonical link

The log-mass is \\y\log(\mu/(1-\mu)) + \log(1-\mu)\\, linear in \\y\\
once \\\eta = \log(\mu/(1-\mu))\\ is the parameter. The logit is
therefore the canonical link, and on its scale \$\$\dfrac{\partial
\ell}{\partial \eta} = y - \mu, \qquad \dfrac{\partial^2 \ell}{\partial
\eta^2} = -\mu(1-\mu).\$\$ The second carries no data, so the observed
and the expected information coincide. Under a probit or a complementary
log-log link they do not, and the choice between Fisher scoring and
Newton's method becomes a real one.

## Derivatives

On the parameter scale, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{y - \mu}{\mu(1-\mu)}, \qquad \dfrac{\partial^2 \ell}{\partial
\mu^2} = -\dfrac{y}{\mu^2} - \dfrac{1-y}{(1-\mu)^2},\$\$ with
expectations \\0\\ and \\-1/(\mu(1-\mu))\\. The information is the
reciprocal of the variance, so it is smallest at \\\mu = 1/2\\ and grows
without bound towards either endpoint: a near-certain outcome is very
informative about its own probability.

## Estimation

The maximum likelihood estimate is the sample proportion, \\\hat\mu =
\bar y\\, in closed form. A sample of all zeros or all ones drives the
estimate to the boundary, where the logit is infinite; that is
separation, and it is a property of the data.

## Notation

\\\ell\\ is the log-mass of one observation and \\\mu \in (0,1)\\ the
success probability. \\\eta = \log(\mu/(1-\mu))\\ is the log odds, the
parameter on the link scale. The **canonical** link of an exponential
family is the one that makes the log-mass linear in the sufficient
statistic.

## References

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). *Univariate Discrete
Distributions*, 3rd edition, Chapter 3. Wiley, Hoboken.

## See also

[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
of which this is the one-trial case;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
and
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for grouped binary data that is overdispersed;
[`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
[`linkfunctions7::probit_link()`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.html)
and
[`linkfunctions7::cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.html)
for the usual links;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the probability;
[BernoulliDistrib](https://statmodels7.github.io/distributions7/reference/BernoulliDistrib.md)
for the class.

## Examples

``` r
d <- bernoulli_distrib()
d
#> Distribution: Bernoulli
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (probability)        | Link: logit      | Domain: (0, 1)

# The two masses, summing to one.
distrib_pdf(d, c(0, 1), list(mu = 0.3))
#> [1] 0.7 0.3

# Mean and variance are tied, the variance being maximal at 1/2.
vapply(c(0.1, 0.5, 0.9), function(p) {
  th <- list(mu = p)
  c(mean = mean(d, th), var = variance(d, th), skew = skewness(d, th))
}, numeric(3))
#>          [,1] [,2]      [,3]
#> mean 0.100000 0.50  0.900000
#> var  0.090000 0.25  0.090000
#> skew 2.666667 0.00 -2.666667

# The estimate is the sample proportion, in closed form.
set.seed(4)
z <- distrib_rng(d, 1000, list(mu = 0.3))
c(fitted = unname(coef(fit_distrib(d, z))), proportion = mean(z))
#>     fitted proportion 
#>      0.287      0.287 

# On the canonical logit link the observed and expected information agree;
# under a probit they do not.
th <- list(mu = 0.3)
probit <- bernoulli_distrib(link_mu = linkfunctions7::probit_link())
rbind(logit = c(obs = distrib_hessian(d, 1, th, scale = "link")$mu_mu,
                exp = distrib_expected_hessian(d, 1, th,
                                               scale = "link")$mu_mu),
      probit = c(obs = distrib_hessian(probit, 1, th, scale = "link")$mu_mu,
                 exp = distrib_expected_hessian(probit, 1, th,
                                                scale = "link")$mu_mu))
#>               obs        exp
#> logit  -0.2100000 -0.2100000
#> probit -0.7354566 -0.5756674
```
