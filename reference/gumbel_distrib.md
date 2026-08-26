# Gumbel Distribution

Builds the distribution object for the Gumbel (type I extreme value)
family in the form for **maxima**: a location-scale family on the whole
real line with location \\\mu\\ and scale \\\sigma \> 0\\. The returned
object carries closed-form derivatives of the log-density to fourth
order, in the parameters and in the response, observed and expected, and
closed-form moments, so every generic of the toolkit answers without a
numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. The defaults are the identity
for the location, which ranges over the whole line, and the logarithm
for the scale, which keeps it positive at every predictor.

## Usage

``` r
gumbel_distrib(link_mu = identity_link(), link_sigma = log_link())
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

An S7 object of class `GumbelDistrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"gumbel"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params` `c("mu", "sigma")`,
`n_params` `2`, `params_bounds` the list of \\(-\infty, \infty)\\ and
\\(0, \infty)\\, and `link_params` the two links given here.

## The parametrization

With \\z = (y-\mu)/\sigma\\ the density on the whole real line is
\$\$f(y; \mu, \sigma) = \dfrac{1}{\sigma} \exp\left\\-z -
e^{-z}\right\\,\$\$ the distribution function \\F(q) = \exp\\-e^{-z}\\\\
and the quantile function \\Q(p) = \mu - \sigma\log(-\log p)\\, both
closed form. The mean is \\\mu + \gamma\sigma\\ with \\\gamma\\ the
Euler-Mascheroni constant, the variance \\\pi^2\sigma^2/6\\, the
skewness \\12\sqrt{6}\\\zeta(3)/\pi^3 \approx 1.1395\\ and the excess
kurtosis \\12/5\\. The last two are constants: **this family has a fixed
shape and only its location and spread can be fitted.** The mode is
\\\mu\\, the median \\\mu - \sigma\log\log 2\\, and the mean is above
both.

The law is the limit of the maximum of a sample from a light-tailed
distribution, once centered and scaled, and that is what it is fitted
to. It is max-stable: raising its distribution function to the \\n\\th
power gives the same law with the location shifted by \\\sigma\log n\\.

## Everything rests on w being standard exponential

Write \\w = e^{-z}\\. Under the model \\w\\ is standard exponential
whatever \\\mu\\ and \\\sigma\\ are, and the whole family is written in
it. The score is \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{1 -
w}{\sigma}, \qquad \dfrac{\partial \ell}{\partial \sigma} = \dfrac{z(1 -
w) - 1}{\sigma},\$\$ the observed Hessian is \$\$\ell^{(\mu\mu)} =
-\dfrac{w}{\sigma^2}, \quad \ell^{(\mu\sigma)} = -\dfrac{1 - w +
zw}{\sigma^2}, \quad \ell^{(\sigma\sigma)} = \dfrac{1 - 2z + 2zw - z^2
w}{\sigma^2},\$\$ and every expectation the family needs is a derivative
of \\\Gamma\\ at 2: \\\mathbb{E}\[w\] = 1\\, \\\mathbb{E}\[w\log w\] =
1-\gamma\\, \\\mathbb{E}\[w(\log w)^2\] = (1-\gamma)^2 + \pi^2/6 - 1\\.
The expected Hessian follows in one line, and so do the third and fourth
orders, through \\\mathbb{E}\[z^k w\] = (-1)^k\Gamma^{(k)}(2)\\.

**The location and the scale are not orthogonal.** The mixed entry of
the expected information is \\(1-\gamma)/\sigma^2\\, which is about
\\0.4228/\sigma^2\\ and never zero. In a symmetric location-scale family
such as the Gaussian it vanishes; here the density is skewed and it does
not, so the two estimates are asymptotically correlated.

## Relatives

If \\Y\\ is Gumbel then \\e^{-Y}\\ is Weibull, so
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
is this family on the log scale and reversed, and the two share the
expectations that produce their information matrices. For **minima** the
law is the reflection: fit this family to \\-Y\\, or use
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
with
[`affine_transform(scale = -1)`](https://statmodels7.github.io/distributions7/reference/affine_transform.md).
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
is the other limit law of extremes, for exceedances over a threshold in
place of block maxima.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Neither estimate is
closed form. The fixed shape supplies the method of moments starting
values: the scale is \\s\sqrt{6}/\pi\\ and the location \\\bar y -
\gamma s\sqrt{6}/\pi\\, with \\s\\ the sample standard deviation. The
example below shows them landing beside the estimates.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location and
\\\sigma \> 0\\ the scale. \\z = (y-\mu)/\sigma\\ is the standardized
value and \\w = e^{-z}\\, standard exponential under the model.
\\\gamma\\ is the Euler-Mascheroni constant, \\-\psi(1) \approx
0.5772\\. \\\eta\\ is a parameter on the unconstrained scale of its
link, with \\\theta = g^{-1}(\eta)\\.

## References

Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme
Values*, Chapter 3. Springer, London.

## See also

[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
for the family this becomes under \\e^{-Y}\\;
[`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md)
for exceedances over a threshold;
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
for the reflection that gives minima;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the symmetric location-scale family whose parameters are orthogonal;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[GumbelDistrib](https://statmodels7.github.io/distributions7/reference/GumbelDistrib.md)
for the class.

## Examples

``` r
d <- gumbel_distrib()
d
#> Distribution: Gumbel
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (location)           | Link: identity   | Domain: (-Inf, Inf)
#>   sigma (scale)              | Link: log        | Domain: (0, Inf)

# The density, written out.
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)
all.equal(distrib_pdf(d, y, th), exp(-y - exp(-y)))
#> [1] TRUE

# The mean is shifted by Euler's constant and the shape is fixed: the
# skewness and kurtosis are the same at any parameter setting.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>      mean       var      skew      kurt 
#> 0.5772157 1.6449341 1.1395471 2.4000000 
c(skewness(d, list(mu = 5, sigma = 3)), kurtosis(d, list(mu = 5, sigma = 3)))
#> [1] 1.139547 2.400000

# Max-stable: the cdf to the nth power is the same law shifted by
# sigma log n.
all.equal(distrib_cdf(d, y, th)^10,
          distrib_cdf(d, y, list(mu = log(10), sigma = 1)))
#> [1] TRUE

# exp(-Y) is Weibull with scale exp(-mu) and shape 1/sigma; the densities
# agree once the Jacobian is applied.
set.seed(3)
yy <- distrib_rng(d, 5, list(mu = 0, sigma = 1 / 2))
all.equal(distrib_pdf(d, yy, list(mu = 0, sigma = 1 / 2)),
          distrib_pdf(weibull1_distrib(), exp(-yy),
                      list(mu = 1, sigma = 2)) * exp(-yy))
#> [1] TRUE

# The location and the scale are not orthogonal: the mixed entry of the
# expected information is (1 - gamma)/sigma^2.
c(mixed = distrib_expected_hessian(d, 0, th)$mu_sigma,
  one_minus_gamma = 1 + digamma(1))
#>           mixed one_minus_gamma 
#>       0.4227843       0.4227843 

# Fitting recovers the parameters; the moment estimates start it off.
set.seed(8)
s <- distrib_rng(d, 2000, list(mu = 3, sigma = 2))
sc <- sd(s) * sqrt(6) / pi
rbind(fitted  = coef(fit_distrib(d, s)),
      moments = c(mu = mean(s) + digamma(1) * sc, sigma = sc))
#>               mu    sigma
#> fitted  2.992365 1.984242
#> moments 2.992978 1.970070
```
