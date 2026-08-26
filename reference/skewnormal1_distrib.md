# Skew Normal Distribution Object

Builds a skew normal distribution object in Azzalini's direct
parametrization: a location \\\mu\\, a scale \\\sigma \> 0\\ and a shape
\\\alpha\\ that tilts the Gaussian. With \\z = (y-\mu)/\sigma\\ the
density is \\2\phi(z)\Phi(\alpha z)/\sigma\\, and \\\alpha = 0\\ gives
the Gaussian exactly.

The returned object carries the three link functions and every method of
the family. `mu` and `sigma` are a location and a scale, not the mean
and the standard deviation;
[`mean.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal1Distrib.md)
and
[`variance.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal1Distrib.md)
give those.

## Usage

``` r
skewnormal1_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_alpha = identity_link()
)
```

## Arguments

- link_mu:

  A `linkfunctions7` link object for the location \\\mu\\, which is
  unconstrained. Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link object for the scale \\\sigma\\, which must be strictly
  positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  so that any real linear predictor maps to an admissible scale.

- link_alpha:

  A link object for the shape \\\alpha\\, which is unconstrained.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

## Value

An S7 object of class
[SkewNormal1Distrib](https://statmodels7.github.io/distributions7/reference/SkewNormal1Distrib.md),
inheriting from `continuous_distrib`. Its `params` are
`c("mu", "sigma", "alpha")`, its `bounds` `c(-Inf, Inf)`, and its
`link_params` the three links given here.

## Density and distribution function

\$\$f(y; \mu, \sigma, \alpha) = \dfrac{2}{\sigma}\\\phi(z)\\\Phi(\alpha
z), \qquad z = (y-\mu)/\sigma.\$\$ The factor \\2\Phi(\alpha z)\\
exceeds one where \\\alpha z \> 0\\ and falls below it where \\\alpha z
\< 0\\, so a positive \\\alpha\\ moves mass to the right. The
distribution function is Azzalini's identity \\F(q) = \Phi(z) - 2T(z,
\alpha)\\ with \\T\\ Owen's T, one bounded quadrature per observation
through
[`numericals7::owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.html).
The quantile function has no closed form and comes from the base class
by root finding.

## Derivatives

The score and the observed Hessian are closed form, written in the
inverse Mills ratio \\R(t) = \phi(t)/\Phi(t)\\ at \\t = \alpha z\\, and
so are the third and fourth orders, in compiled kernels. One identity
does all of it: \\R' = -R(t+R)\\, so every derivative of \\\log\Phi(t)\\
is a polynomial in \\t\\ and \\R\\.

The **expected** information has no elementary form, so none is
registered and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
approximates it by the strategy named in its `approx` argument. The
expected third and fourth orders share that obstruction.

## Singularity at symmetry

At \\\alpha = 0\\ the expected information has rank 2 and not 3. The
score shows why: the shape component is \\z\sqrt{2/\pi}\\ and the
location component \\z/\sigma\\, so the two are exactly proportional and
no data can separate them. Measured, the smallest eigenvalue is
\\-5.6\times10^{-17}\\ against a largest of 2 at \\\alpha = 0\\, and it
grows like \\\alpha^4\\ thereafter.

The consequence for use is that a fit whose true shape is near zero
identifies \\\alpha\\ weakly, and that a variance matrix computed at
exactly zero is not invertible. Azzalini and Capitanio's centered
parametrization removes the singularity; it is
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md),
a separate object.

## Moments, and the bound on the skewness

With \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and \\b = \sqrt{2/\pi}\\,
\$\$E\[Y\] = \mu + \sigma b \delta, \qquad \mathrm{Var}(Y) =
\sigma^2(1 - b^2\delta^2).\$\$ All four moments are closed form. The
skewness is bounded: \\\delta\\ saturates at one as
\\\|\alpha\|\to\infty\\, so \\\gamma_1\\ cannot leave \\(-0.9953,
0.9953)\\ whatever the shape is. Measured, it reaches 0.9556 at \\\alpha
= 10\\ and 0.99527 at \\\alpha = 10^4\\, and does not move after that. A
sample skewer than this needs
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

## Parameter domains

- \\\mu \in (-\infty, \infty)\\

- \\\sigma \in (0, \infty)\\

- \\\alpha \in (-\infty, \infty)\\

## References

Azzalini, A. (1985). A class of distributions which includes the normal
ones. *Scandinavian Journal of Statistics* 12, 171-178.

Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
Families*. Cambridge University Press.

## See also

[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
for the centered parametrization,
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
for the four-parameter family that reaches a larger skewness,
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the case \\\alpha = 0\\, and
[SkewNormal1Distrib](https://statmodels7.github.io/distributions7/reference/SkewNormal1Distrib.md)
for the class and its method list.

## Examples

``` r
d <- skewnormal1_distrib()
d@params
#> [1] "mu"    "sigma" "alpha"
th <- list(mu = 0, sigma = 1, alpha = 3)

# The location is not the mean, and the scale is not the standard deviation.
rbind(parameter = c(th$mu, th$sigma),
      moment = c(mean(d, th), sqrt(variance(d, th))))
#>                [,1]      [,2]
#> parameter 0.0000000 1.0000000
#> moment    0.7569398 0.6534847

# Shape zero is the Gaussian, in the density and in the distribution function.
y <- c(-1, 0, 1)
th0 <- list(mu = 0, sigma = 1, alpha = 0)
c(density = max(abs(distrib_pdf(d, y, th0) - dnorm(y))),
  cdf = max(abs(distrib_cdf(d, y, th0) - pnorm(y))))
#>      density          cdf 
#> 5.551115e-17 0.000000e+00 

# The skewness the family can reach is bounded, and saturates early.
vapply(c(1, 3, 10, 50, 1e4),
       function(a) skewness(d, list(mu = 0, sigma = 1, alpha = a)), 0)
#> [1] 0.1369488 0.6670236 0.9555571 0.9936306 0.9952717

# A fit recovers all three parameters.
set.seed(5)
x <- distrib_rng(d, 3000, list(mu = 2, sigma = 1.5, alpha = 4))
coef(fit_distrib(d, x))
#>       mu    sigma    alpha 
#> 2.010681 1.505278 3.963568 
```
