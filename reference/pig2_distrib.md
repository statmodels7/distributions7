# Poisson-Inverse Gaussian Distribution Object, Orthogonal Parametrization

Builds a Poisson-inverse Gaussian distribution in the parametrization
whose two parameters are orthogonal, which is gamlss's `PIG2`. The mean
stays \\\mu\\; the second parameter \\\alpha\\ is exactly the argument
at which the mass function's Bessel function is evaluated, related to
the dispersion of
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
by \\\alpha = \sqrt{1 + 2\sigma\mu}/\sigma\\.

Orthogonal means the expected information is diagonal. Measured at \\\mu
= 3\\, its mixed entry summed over the support is \\-8.8\times10^{-15}\\
here against 7.39 in the mean-dispersion parametrization.

## Usage

``` r
pig2_distrib(link_mu = log_link(), link_alpha = log_link())
```

## Arguments

- link_mu:

  A `linkfunctions7` link object for the mean \\\mu\\, which must be
  strictly positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_alpha:

  A link object for \\\alpha\\, also strictly positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class
[Pig2Distrib](https://statmodels7.github.io/distributions7/reference/Pig2Distrib.md),
inheriting from `discrete_distrib`. Its `params` are `c("mu", "alpha")`,
its `bounds` `c(0, Inf)`, and its `link_params` the two links given
here.

## What orthogonality buys, and what it costs

The maximum likelihood estimates of \\\mu\\ and \\\alpha\\ are
asymptotically independent, so a Fisher scoring step in one parameter
does not disturb the other and the two can be modeled with separate
linear predictors without their estimates fighting. The cost is that
\\\alpha\\ has no moment reading of its own: a reader who wants to
interpret the overdispersion directly wants
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)'s
\\\sigma\\, which enters the variance as \\\mu + \sigma\mu^2\\.

The two are one law and the map between them is closed both ways:
\\\sigma = (\mu + \sqrt{\mu^2 + \alpha^2})/\alpha^2\\ through
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md),
and \\\alpha = \sqrt{1 + 2\sigma\mu}/\sigma\\ back. gamlss states the
first as \\\alpha = 1/(\sqrt{\mu^2+\sigma_2^2} - \mu)\\ in its own
\\\sigma_2\\, which coincides with this \\\alpha\\.

## How the derivatives are computed

By the same compiled kernel as
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
with \\\alpha\\ as a variable of its own, so the Bessel argument needs
no chain rule at all. Every partial to fourth order is exact and comes
out of one pass; see
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md).
The expected information has no closed form and goes through
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

## Parameter domains

- \\\mu \in (0, \infty)\\

- \\\alpha \in (0, \infty)\\

A large \\\alpha\\ is a small dispersion, and the family tends to the
Poisson: measured at \\\alpha = 10^4\\ the mass agrees with `dpois` over
the head of the support.

## The distribution

\$\$P(Y=y) =
\sqrt{\frac{2\alpha}{\pi}}\\\frac{\mu^{y}e^{1/\sigma}}{(\alpha\sigma)^{y}\\y!}\\K\_{y-1/2}(\alpha),
\qquad \sigma = \frac{1}{\sqrt{\mu^{2}+\alpha^{2}} - \mu}\$\$ on \\y \in
\\0, 1, \dots\\\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \mu +
\sigma\mu^{2}\$\$

## References

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics* 54(3),
507–554.

Heller, G. Z., Couturier, D.-L., and Heritier, S. R. (2019). Beyond mean
modelling: bias due to misspecification of dispersion in Poisson-inverse
Gaussian regression. *Biometrical Journal* 61(2), 333–342.

## See also

[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the mean-dispersion parametrization,
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
for the map,
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the other overdispersed count family, and
[Pig2Distrib](https://statmodels7.github.io/distributions7/reference/Pig2Distrib.md)
for the class.

## Examples

``` r
d <- pig2_distrib()
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
th <- list(mu = 3, alpha = al)

distrib_pdf(d, 0:5, th)
#> [1] 0.17197629 0.21422781 0.17775288 0.12895666 0.08968701 0.06196187
c(mean = mean(d, th), variance = variance(d, th))
#>     mean variance 
#>      3.0     10.2 

# The same law as pig1 at the dispersion alpha implies.
all.equal(distrib_pdf(d, 0:5, th),
          distrib_pdf(pig1_distrib(), 0:5, list(mu = 3, sigma = 0.8)))
#> [1] TRUE

# The property the parametrization exists for.
c(pig2 = sum(distrib_expected_hessian(d, 0:200, th,
                                      approx = "bartlett")$mu_alpha),
  pig1 = sum(distrib_expected_hessian(pig1_distrib(), 0:200,
                                      list(mu = 3, sigma = 0.8),
                                      approx = "bartlett")$mu_sigma))
#>          pig2          pig1 
#> -1.456596e-14  7.392208e+00 

# A large alpha is a small dispersion, and the family tends to the Poisson.
rbind(pig2 = distrib_pdf(d, 0:5, list(mu = 3, alpha = 1e4)),
      poisson = dpois(0:5, 3))
#>               [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> pig2    0.04980948 0.1493836 0.2240306 0.2240082 0.1680061 0.1008138
#> poisson 0.04978707 0.1493612 0.2240418 0.2240418 0.1680314 0.1008188

# A fit recovers both parameters.
set.seed(65)
x <- distrib_rng(d, 4000, th)
coef(fit_distrib(d, x))
#>       mu    alpha 
#> 3.041000 2.983942 
```
