# Gamma Distribution, Mean and Dispersion

Builds the distribution object for the gamma family parametrized by its
mean \\\mu \> 0\\ and a dispersion \\\phi \> 0\\, so that
\\\operatorname{Var}(Y) = \phi\mu^2\\. This is the parametrization a
generalized linear model uses. The returned object carries closed-form
derivatives of the log-density to fourth order, in the parameters and in
the response, and closed-form moments, so every generic of the toolkit
answers without a numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both parameters being positive.

## Usage

``` r
gamma1_distrib(link_mu = log_link(), link_phi = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  mean positive.

- link_phi:

  A `link` object from `linkfunctions7` for the dispersion \\\phi\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `Gamma1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"gamma1"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "phi")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## The parametrization

The density on \\y \in (0, \infty)\\ is \$\$f(y; \mu, \phi) =
\dfrac{y^{1/\phi - 1}\\e^{-y/(\phi\mu)}}
{(\phi\mu)^{1/\phi}\\\Gamma(1/\phi)},\$\$ the gamma at shape \\a =
1/\phi\\ and rate \\b = 1/(\phi\mu)\\. The mean is \\\mu\\, the variance
\\\phi\mu^2\\, the skewness \\2\sqrt{\phi}\\ and the excess kurtosis
\\6\phi\\. The coefficient of variation is \\\sqrt{\phi}\\ at every
mean, so \\\phi\\ measures relative rather than absolute spread.

The variance function is \\V(\mu) = \mu^2\\ and \\\phi\\ is the
dispersion multiplying it, which makes the score in \\\mu\\ the
generalized linear model score \\(y-\mu)/(\phi\mu^2)\\. At \\\phi = 1\\
the shape is 1 and the family is the exponential,
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md);
at \\\phi \to 0\\ it tends to a Gaussian with variance \\\phi\mu^2\\.

This is the same law as
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
which carries the mean and the *variance*, the two being related by
\\\sigma^2 = \phi\mu^2\\. They are separate families because the second
parameter is a different quantity in each, with its own interpretation,
standard error and interval.

## Derivatives

Writing \\s = 1/\phi\\ and \\z = y/\mu\\, the score is
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y-\mu}{\phi\mu^2},
\qquad \dfrac{\partial \ell}{\partial \phi} = -s^2\left\\\log s + 1 -
\psi(s) + \log z - z\right\\,\$\$ and the expected Hessian is
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{1}{\phi\mu^2},
\quad \mathbb{E}\left\[\ell^{(\mu\phi)}\right\] = 0, \quad
\mathbb{E}\left\[\ell^{(\phi\phi)}\right\] = s^4\left\\\dfrac{1}{s} -
\psi'(s)\right\\.\$\$ The zero off-diagonal makes the mean and the
dispersion orthogonal, so the mean equation can be fitted with the
dispersion held at any value without biasing it.

Every derivative in \\\phi\\ is the corresponding derivative in \\s\\
carried across by the one-variable chain rule, with \\s' = -s^2\\, \\s''
= 2s^3\\, \\s''' = -6s^4\\ and \\s'''' = 24s^5\\, so each polygamma
function is evaluated once. Two quantities in those derivatives cancel
as \\s\\ grows, \\\log s - \psi(s)\\ and its polygamma analogues; each
is computed as a polygamma minus its own leading asymptote, so the
digits survive at the large shape a nearly Gaussian gamma reaches.

Third and fourth orders are closed form as well, observed and expected,
in
[`distrib_deriv3.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma1Distrib.md)
and
[`distrib_deriv4.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gamma1Distrib.md),
as are the derivatives in the response,
[`distrib_grad_y.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.Gamma1Distrib.md)
and
[`distrib_hess_y.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.Gamma1Distrib.md).
The derivatives of the *distribution* function in the parameters have no
elementary form here, the derivative of an incomplete gamma in its shape
being hypergeometric, and are taken by finite difference on the analytic
cdf.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. The mean has the
closed-form estimate \\\hat\mu = \bar y\\. The shape \\\hat s\\ then
solves \$\$\log \hat s - \psi(\hat s) = -\dfrac{1}{n}\sum_i
\log(y_i/\bar y),\$\$ which has no closed form and is reached
numerically, and \\\hat\phi = 1/\hat s\\. The method of moments starting
value is the squared coefficient of variation, and the example below
shows it landing beside the estimate.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the mean
and \\\phi \> 0\\ the dispersion. \\\psi\\ and \\\psi'\\ are the digamma
and trigamma functions. \\\eta\\ is a parameter on the unconstrained
scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
edition, Chapter 8. Chapman and Hall, London.

## See also

[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
for the same law in the mean and the variance;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the case \\\phi = 1\\;
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md)
and
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
for relatives;
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
and
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for other positive families with a multiplicative variance function;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Gamma1Distrib](https://statmodels7.github.io/distributions7/reference/Gamma1Distrib.md)
for the class.

## Examples

``` r
d <- gamma1_distrib()
d
#> Distribution: Gamma1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu  (mean)               | Link: log        | Domain: (0, Inf)
#>   phi (dispersion)         | Link: log        | Domain: (0, Inf)

# The density is stats::dgamma at shape 1/phi and rate 1/(phi mu).
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)
all.equal(distrib_pdf(d, y, th), dgamma(y, shape = 2, rate = 1 / 1.5))
#> [1] TRUE

# Moments in closed form: skewness 2 sqrt(phi), excess kurtosis 6 phi.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>     mean      var     skew     kurt 
#> 3.000000 4.500000 1.414214 3.000000 
c(2 * sqrt(0.5), 6 * 0.5)
#> [1] 1.414214 3.000000

# Fitting recovers the parameters, the mean exactly at the sample mean.
set.seed(4)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted  = coef(fit),
      moments = c(mu = mean(z), phi = var(z) / mean(z)^2))
#>               mu       phi
#> fitted  2.922086 0.4959732
#> moments 2.922086 0.5015235

# At phi = 1 the family is the exponential.
all.equal(distrib_pdf(d, y, list(mu = 3, phi = 1)),
          distrib_pdf(exponential_distrib(), y, list(mu = 3)))
#> [1] TRUE
```
