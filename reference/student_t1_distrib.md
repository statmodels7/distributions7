# Student t Distribution, Location, Scale and Degrees of Freedom

Builds the distribution object for the location-scale Student t family,
parametrized by a location \\\mu\\, a scale \\\sigma \> 0\\ and degrees
of freedom \\\nu \> 0\\. The returned object carries closed-form
derivatives of the log-density to fourth order in the parameters, closed
first and second derivatives in the response, and a closed expected
Hessian; the expected third and fourth orders are the only quantities
that go through a numerical route.

The family is the standard heavy-tailed alternative to a Gaussian. Its
location score redescends, so a gross outlier contributes almost nothing
to the estimating equation, and \\\nu\\ is estimated from the data
rather than set: a large fitted \\\nu\\ reports that the sample looks
Gaussian.

The three arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on.

## Usage

``` r
student_t1_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_nu = log_link()
)
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

- link_nu:

  A `link` object from `linkfunctions7` for the degrees of freedom
  \\\nu\\. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `StudentT1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"student t1"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params`
`c("mu", "sigma", "nu")`, `n_params` `3`, `params_bounds` the domains
\\(-\infty, \infty)\\, \\(0, \infty)\\ and \\(0, \infty)\\, and
`link_params` the three links given here.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \sigma, \nu)
=
\dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\\\Gamma\left(\dfrac{\nu}{2}\right)}
\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}},\$\$
the distribution function \\F(q) = T\_\nu\\(q-\mu)/\sigma\\\\ and the
quantile function \\Q(p) = \mu + \sigma T\_\nu^{-1}(p)\\, with
\\T\_\nu\\ the standard Student t distribution function.

\\\sigma\\ is the **scale** and not the standard deviation, which is
\\\sigma\sqrt{\nu/(\nu-2)}\\ and exists only above \\\nu = 2\\. Keeping
the two apart lets the family be fitted where the second moment does not
exist;
[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md)
is the parametrization by the standard deviation, for a reader who wants
one.

Two shapes are named families: \\\nu = 1\\ is the Cauchy, and \\\nu \to
\infty\\ the Gaussian with standard deviation \\\sigma\\. The approach
to the limit is \\O(1/\nu)\\, so at \\\nu = 10^5\\ the density already
agrees with a Gaussian's to seven figures.

## Derivatives

Write \\r = y - \mu\\ and \\D = \nu\sigma^2 + r^2\\. The score is
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{(\nu+1)r}{D}, \qquad
\dfrac{\partial \ell}{\partial \sigma} = \dfrac{\nu\left(r^2 -
\sigma^2\right)}{\sigma D},\$\$ \$\$\dfrac{\partial \ell}{\partial \nu}
= \dfrac{1}{2}\left\[ -\dfrac{1}{\nu} -
\psi\\\left(\dfrac{\nu}{2}\right) +
\psi\\\left(\dfrac{\nu+1}{2}\right) + \dfrac{(\nu+1)r^2}{\nu D} -
\log\\\left(1 + \dfrac{r^2}{\nu\sigma^2}\right)\right\],\$\$ with
\\\psi\\ the digamma function. The observed Hessian is in
[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md),
and its expectation is \$\$\mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu^2}\right\] = -\dfrac{\nu+1}{\sigma^2(\nu+3)}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{2\nu}{\sigma^2(\nu+3)},\$\$
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \nu^2}\right\] =
\dfrac{1}{4}\left\[\psi_1\\\left(\dfrac{\nu+1}{2}\right) -
\psi_1\\\left(\dfrac{\nu}{2}\right)\right\] +
\dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma \\ \partial
\nu}\right\] = \dfrac{2}{\sigma(\nu+1)(\nu+3)},\$\$ with the two entries
containing \\\mu\\ exactly zero. The location is therefore orthogonal to
the other two, whose estimates are asymptotically correlated with each
other.

## Redescent, and where it costs

The location score \\(\nu+1)r/D\\ rises to
\\(\nu+1)/(2\sigma\sqrt{\nu})\\ at \\\|r\| = \sigma\sqrt{\nu}\\ and
falls back towards zero, so an observation far from the location is
downweighted automatically. The price is that the curvature in \\\mu\\
turns **positive** past the same point, so the observed information is
indefinite at an outlying observation. Fisher scoring is unaffected, the
expected information being negative definite everywhere, which is why it
is the default in
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md).

## Moments

Each exists only above its own threshold, and the family reports `NaN`
or `Inf` where it does not: the mean is \\\mu\\ for \\\nu \> 1\\, the
variance \\\sigma^2\nu/(\nu-2)\\ for \\\nu \> 2\\, the skewness 0 for
\\\nu \> 3\\ and the excess kurtosis \\6/(\nu-4)\\ for \\\nu \> 4\\. At
\\\nu = 1.5\\, for instance, the mean is \\\mu\\ and the variance is
`Inf`.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale by Fisher scoring. No
estimate is closed form, the three estimating equations being coupled
through \\D\\. On Gaussian data \\\hat\nu\\ runs towards its upper
boundary, which is the correct answer and is reported as a value near
the clamp of its log link rather than as `Inf`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location,
\\\sigma \> 0\\ the scale and \\\nu \> 0\\ the degrees of freedom.
\\\psi\\ and \\\psi_1\\ are the digamma and trigamma functions,
[`digamma()`](https://rdrr.io/r/base/Special.html) and
[`trigamma()`](https://rdrr.io/r/base/Special.html) in R. \\\eta\\ is a
parameter on the unconstrained scale of its link, with \\\theta =
g^{-1}(\eta)\\.

## References

Lange, K. L., Little, R. J. A. and Taylor, J. M. G. (1989). Robust
statistical modeling using the t distribution. *Journal of the American
Statistical Association*, **84**(408), 881-896.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995). *Continuous
Univariate Distributions*, Volume 2, 2nd edition, Chapter 28. Wiley, New
York.

## See also

[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md)
for the same law in the standard deviation;
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
for \\\nu = 1\\ and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the limit;
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
to add a shape parameter;
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
and
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
for other robust families;
[`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
for the multivariate version;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[StudentT1Distrib](https://statmodels7.github.io/distributions7/reference/StudentT1Distrib.md)
for the class.

## Examples

``` r
d <- student_t1_distrib()
d
#> Distribution: Student T1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (location)           | Link: identity   | Domain: (-Inf, Inf)
#>   sigma (scale)              | Link: log        | Domain: (0, Inf)
#>   nu    (shape)              | Link: log        | Domain: (0, Inf)

# The density is stats::dt at the standardized value, over sigma.
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)
all.equal(distrib_pdf(d, y, th), dt((y - 0.4) / 1.2, df = 5) / 1.2)
#> [1] TRUE

# The scale is not the standard deviation; the moments carry thresholds.
c(scale = th$sigma, mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#> scale  mean   var  skew  kurt 
#>   1.2   0.4   2.4   0.0   6.0 

# Below its own threshold each moment reports that it does not exist.
t(vapply(c(0.8, 1.5, 2.5, 3.5),
         function(v) {
           p <- list(mu = 0.4, sigma = 1.2, nu = v)
           c(nu = v, mean = mean(d, p), var = variance(d, p),
             kurt = kurtosis(d, p))
         }, numeric(4)))
#>       nu mean  var kurt
#> [1,] 0.8  NaN  NaN  NaN
#> [2,] 1.5  0.4  Inf  NaN
#> [3,] 2.5  0.4 7.20  Inf
#> [4,] 3.5  0.4 3.36  Inf

# One degree of freedom is the Cauchy; a large one is the Gaussian.
all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 1)),
          dcauchy(y, 0.4, 1.2))
#> [1] TRUE
max(abs(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 1e5)) -
        dnorm(y, 0.4, 1.2)))
#> [1] 9.603903e-07

# Fitting recovers all three parameters.
set.seed(4)
z <- distrib_rng(d, 3000, list(mu = 1, sigma = 2, nu = 4))
coef(fit_distrib(d, z))
#>       mu    sigma       nu 
#> 1.021176 1.967036 3.998383 
```
