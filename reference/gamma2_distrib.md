# Gamma Distribution, Mean and Variance

Builds the distribution object for the gamma family parametrized by its
mean \\\mu \> 0\\ and its variance \\\sigma^2 \> 0\\, so that the shape
is \\\alpha = \mu^2/\sigma^2\\ and the rate \\\lambda = \mu/\sigma^2\\.
The returned object carries closed-form derivatives of the log-density
to fourth order, in the parameters and in the response, and closed-form
moments, so every generic of the toolkit answers without a numerical
fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both parameters being positive.

## Usage

``` r
gamma2_distrib(link_mu = log_link(), link_sigma2 = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  mean positive.

- link_sigma2:

  A `link` object from `linkfunctions7` for the variance \\\sigma^2\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `Gamma2Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"gamma2"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "sigma2")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## The parametrization

The density on \\y \in (0, \infty)\\ is \$\$f(y; \mu, \sigma^2) =
\dfrac{\lambda^{\alpha}}{\Gamma(\alpha)}\\ y^{\alpha-1} e^{-\lambda y},
\qquad \alpha = \dfrac{\mu^2}{\sigma^2}, \quad \lambda =
\dfrac{\mu}{\sigma^2},\$\$ the distribution function \\F(q) =
\gamma(\alpha, \lambda q)/\Gamma(\alpha)\\ with \\\gamma\\ the lower
incomplete gamma function, and the quantile function its numerical
inverse. The mean is \\\mu\\, the variance \\\sigma^2\\, the skewness
\\2\sqrt{\sigma^2}/\mu\\ and the excess kurtosis \\6\sigma^2/\mu^2\\.

This is the same law as
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md),
which carries the mean and a *dispersion*, the two being related by
\\\sigma^2 = \phi\mu^2\\. They are separate families because the second
parameter is a different quantity in each, with its own interpretation,
standard error and interval. The choice between them is not only
cosmetic: see the next section.

## Derivatives, and orthogonality

Writing \\\psi\\ for the digamma function and \\\psi_1\\ for the
trigamma, the score is \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{-2\mu\psi(\alpha) + 2\mu\log\lambda + \mu + 2\mu\log y - y}
{\sigma^2}, \qquad \dfrac{\partial \ell}{\partial \sigma^2} =
-\dfrac{\mu\left\\-\mu\psi(\alpha) + \mu + \mu(\log\lambda + \log y) -
y\right\\}{(\sigma^2)^2},\$\$ and the expected Hessian is
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = \dfrac{3\sigma^2 -
4\mu^2\psi_1(\alpha)}{(\sigma^2)^2}, \qquad
\mathbb{E}\left\[\ell^{(\sigma^2\sigma^2)}\right\] =
-\dfrac{\mu^2\left\\\mu^2\psi_1(\alpha) - \sigma^2\right\\}
{(\sigma^2)^4},\$\$ \$\$\mathbb{E}\left\[\ell^{(\mu\sigma^2)}\right\] =
\dfrac{2\mu\left\\\mu^2\psi_1(\alpha) - \sigma^2\right\\}
{(\sigma^2)^3}.\$\$

**The mixed entry does not vanish.** The mean and the variance are not
orthogonal here, so their estimates are asymptotically correlated and a
mean equation fitted with the variance held at a wrong value is biased.
In
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md),
where the second parameter is the dispersion \\\phi = \sigma^2/\mu^2\\,
the same entry is exactly zero. That is the reason a generalized linear
model uses the dispersion.

Third and fourth orders are closed form as well, observed and expected,
in
[`distrib_deriv3.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma2Distrib.md)
and
[`distrib_deriv4.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gamma2Distrib.md),
as are the derivatives in the response,
[`distrib_grad_y.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gamma2Distrib.md)
and
[`distrib_hess_y.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gamma2Distrib.md).
The derivatives of the *distribution* function in the parameters have no
elementary form, the derivative of an incomplete gamma in its shape
being hypergeometric, and are taken by finite difference on the analytic
cdf.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. The mean has the
closed-form estimate \\\hat\mu = \bar y\\; the variance has none and is
reached numerically, landing near but not at the sample variance. The
example below shows both.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the mean
and \\\sigma^2 \> 0\\ the variance. \\\alpha\\ and \\\lambda\\ are the
implied shape and rate. \\\psi\\ is the digamma function and \\\psi_m\\
its \\m\\th derivative. \\\eta\\ is a parameter on the unconstrained
scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994). *Continuous
Univariate Distributions*, Volume 1, 2nd edition, Chapter 17. Wiley, New
York.

## See also

[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
for the same law in the mean and the dispersion, which is the orthogonal
parametrization;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the case \\\sigma^2 = \mu^2\\;
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md)
and
[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md)
for other positive families written in the mean and the variance;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Gamma2Distrib](https://statmodels7.github.io/distributions7/reference/Gamma2Distrib.md)
for the class.

## Examples

``` r
d <- gamma2_distrib()
d
#> Distribution: Gamma2
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu     (mean)               | Link: log        | Domain: (0, Inf)
#>   sigma2 (variance)           | Link: log        | Domain: (0, Inf)

# The density is stats::dgamma at shape mu^2/sigma2 and rate mu/sigma2.
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)
all.equal(distrib_pdf(d, y, th), dgamma(y, shape = 9 / 2, rate = 3 / 2))
#> [1] TRUE

# Moments: skewness 2 sqrt(sigma2)/mu, excess kurtosis 6 sigma2/mu^2.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>     mean      var     skew     kurt 
#> 3.000000 2.000000 0.942809 1.333333 
c(2 * sqrt(2) / 3, 6 * 2 / 9)
#> [1] 0.942809 1.333333

# Fitting recovers the parameters. The mean is the sample mean exactly;
# the variance is close to the sample variance without equalling it.
set.seed(5)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted  = coef(fit),
      moments = c(mu = mean(z), sigma2 = var(z)))
#>               mu   sigma2
#> fitted  3.025125 2.052645
#> moments 3.025125 2.056397

# The mean and the variance are correlated here, so the mixed entry of the
# expected Hessian is non-zero; in gamma1 it is exactly zero.
c(gamma2 = distrib_expected_hessian(d, 0, th)$mu_sigma2,
  gamma1 = distrib_expected_hessian(gamma1_distrib(), 0,
                                    list(mu = 3, phi = 2 / 9))$mu_phi)
#>    gamma2    gamma1 
#> 0.1788944 0.0000000 
```
