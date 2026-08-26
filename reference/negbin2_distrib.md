# Negative Binomial Distribution, NB2

Builds the distribution object for the negative binomial family on the
non-negative integers in the NB2 parametrization: the mean \\\mu \> 0\\
and a dispersion \\\theta \> 0\\, so that \\\operatorname{Var}(Y) =
\mu + \mu^2/\theta\\. This is the count model a regression reaches for
when a Poisson is overdispersed. The returned object carries closed-form
derivatives of the log-mass to fourth order and closed-form moments.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both parameters being positive.

## Usage

``` r
negbin2_distrib(link_mu = log_link(), link_theta = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  mean positive.

- link_theta:

  A `link` object from `linkfunctions7` for the dispersion \\\theta\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `NegBin2Distrib`, inheriting from
`discrete_distrib`, with `distrib_name` `"negbin2"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "theta")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## The parametrization

The mass on \\y = 0, 1, 2, \ldots\\ is \$\$P(Y = y; \mu, \theta) =
\dfrac{\Gamma(y+\theta)}{y!\\\Gamma(\theta)}
\left(\dfrac{\theta}{s}\right)^{\theta} \left(\dfrac{\mu}{s}\right)^{y},
\qquad s = \theta + \mu,\$\$ the distribution function is the partial
sum and the quantile function its generalized inverse. The mean is
\\\mu\\, the variance \\\mu + \mu^2/\theta\\, the skewness \\(\theta +
2\mu)/\sqrt{\mu\theta(\theta+\mu)}\\ and the excess kurtosis
\\6/\theta + \theta/\\\mu(\theta+\mu)\\\\.

The law is a Poisson whose rate is gamma distributed with mean \\\mu\\
and shape \\\theta\\. Two limits follow: \\\theta = 1\\ gives the
geometric,
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md),
and \\\theta \to \infty\\ gives the Poisson,
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).
**The dispersion reads the opposite way round from a variance**: a small
\\\theta\\ is heavy overdispersion and a large one is nearly a Poisson.

The numbering follows Cameron and Trivedi: NB2 has a variance quadratic
in the mean, and
[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
has one linear in it.

## Derivatives

With \\\psi\\ the digamma function, the score is \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{\theta}{s}\left(\dfrac{y}{\mu} - 1\right),
\qquad \dfrac{\partial \ell}{\partial \theta} = \psi(y+\theta) -
\psi(\theta) + \log\dfrac{\theta}{s} + \dfrac{\mu - y}{s},\$\$ and the
expected Hessian is \$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] =
-\dfrac{\theta}{\mu s}, \qquad
\mathbb{E}\left\[\ell^{(\mu\theta)}\right\] = 0, \qquad
\mathbb{E}\left\[\ell^{(\theta\theta)}\right\] =
\mathbb{E}\[\psi_1(Y+\theta)\] - \psi_1(\theta) + \dfrac{\mu}{\theta
s}.\$\$ The zero off-diagonal makes the mean and the dispersion
orthogonal.

## Every derivative in theta cancels as theta grows

The family tends to the Poisson as \\\theta\\ grows, so every derivative
in \\\theta\\ vanishes there and is written as a sum of terms that
cancel to leading order. The score and the observed Hessian are computed
in forms that perform those cancellations symbolically, and are reliable
throughout; the direct expression for the score is wrong by one part in
a thousand at \\\theta = 10^6\\ and changes sign at \\10^8\\. Two
quantities are **not** rewritten, and both pages say so:

- the expected information in \\\theta\\, which is measured positive
  from about \\\theta = 10^6\\ and so is unusable there
  ([`distrib_expected_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin2Distrib.md));

- the third and fourth derivatives in \\\theta\\
  ([`distrib_deriv3.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin2Distrib.md),
  [`distrib_deriv4.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBin2Distrib.md)).

A fit reaches that regime routinely, a nearly equidispersed sample
driving \\\theta\\ towards its boundary.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. The mean has the
closed-form estimate \\\hat\mu = \bar y\\; the dispersion solves an
equation in digamma functions and is reached numerically. The method of
moments supplies the starting value \\\bar y^2/(s^2 - \bar y)\\, with
\\s^2\\ the sample variance.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \> 0\\ the mean and
\\\theta \> 0\\ the dispersion. \\\psi\\ is the digamma function and
\\\psi_m\\ its \\m\\th derivative. \\\eta\\ is a parameter on the
unconstrained scale of its link, with \\\theta_j = g^{-1}(\eta_j)\\.

## References

Cameron, A. C. and Trivedi, P. K. (2013). *Regression Analysis of Count
Data*, 2nd edition, Chapter 3. Cambridge University Press, Cambridge.

## See also

[`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md)
for the NB1 parametrization, with a variance linear in the mean;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
for the equidispersed limit and
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
for \\\theta = 1\\;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the Poisson mixed over an inverse Gaussian instead of a gamma;
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for counts with excess zeros;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[NegBin2Distrib](https://statmodels7.github.io/distributions7/reference/NegBin2Distrib.md)
for the class.

## Examples

``` r
d <- negbin2_distrib()
d
#> Distribution: Negbin2
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (mean)               | Link: log        | Domain: (0, Inf)
#>   theta (dispersion)         | Link: log        | Domain: (0, Inf)

# The mass is stats::dnbinom at size = theta, mu = mu.
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 2)
all.equal(distrib_pdf(d, y, th), dnbinom(y, size = 2, mu = 4))
#> [1] TRUE

# Moments: the variance exceeds the mean, by mu^2/theta.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>      mean       var      skew      kurt 
#>  4.000000 12.000000  1.443376  3.083333 
c(4 + 4^2 / 2, 6 / 2 + 2 / (4 * 6))
#> [1] 12.000000  3.083333

# theta = 1 is the geometric and a large theta is the Poisson.
all.equal(distrib_pdf(d, y, list(mu = 4, theta = 1)),
          distrib_pdf(geometric_distrib(), y, list(mu = 4)))
#> [1] TRUE
rbind(negbin = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e6)),
      poisson = dpois(0:4, 3))
#>               [,1]      [,2]      [,3]      [,4]      [,5]
#> negbin  0.04978729 0.1493614 0.2240417 0.2240415 0.1680311
#> poisson 0.04978707 0.1493612 0.2240418 0.2240418 0.1680314

# Fitting recovers the parameters; the moment estimates start it off.
set.seed(5)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted  = coef(fit),
      moments = c(mu = mean(z),
                  theta = mean(z)^2 / (var(z) - mean(z))))
#>             mu    theta
#> fitted  3.8945 2.092325
#> moments 3.8945 1.936378
```
