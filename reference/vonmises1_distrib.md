# von Mises Distribution, Mean Direction and Concentration

Builds the distribution object for the von Mises family, the natural
distribution for an angle, parametrized by a mean direction \\\mu\\ and
a concentration \\\kappa \> 0\\. The returned object carries closed-form
derivatives of the log-density to fourth order in the parameters and in
the response, a closed-form expected information, and a **distribution
function from a Bessel series** in place of the quadrature the base
class would use.

## Usage

``` r
vonmises1_distrib(
  link_mu = bounded_link(lwr = -pi, upr = pi),
  link_kappa = log_link()
)
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean direction \\\mu\\.
  Defaults to `linkfunctions7::bounded_link(lwr = -pi, upr = pi)`, which
  maps the free scale onto \\(-\pi, \pi)\\. See the note on the chart
  below.

- link_kappa:

  A `link` object from `linkfunctions7` for the concentration
  \\\kappa\\. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `VonMises1Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"von mises1"`, `dimension`
`"univariate"`, `bounds` `c(-pi, pi)`, `params` `c("mu", "kappa")`,
`n_params` `2`, `params_bounds` the domains \\(-\pi, \pi)\\ and \\(0,
\infty)\\, and `link_params` the two links given here.

## The parametrization

The density on \\y \in \[-\pi, \pi)\\ is \$\$f(y; \mu, \kappa) =
\dfrac{e^{\kappa\cos(y-\mu)}}{2\pi I_0(\kappa)},\$\$ with \\I_0\\ the
modified Bessel function of the first kind. This is the first family in
the package whose support has the topology of a circle: the two ends of
the interval are the same point, so the density need not vanish at
either.
[`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
is the same law parametrized by the mean resultant length in place of
the concentration.

The normalizing constant goes through
[`numericals7::log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.html).
R's exponentially scaled `besselI` underflows to an exact zero between
\\\kappa = 10^5\\ and \\10^6\\, where the logarithm returns `-Inf`;
measured at \\\kappa = 10^6\\ the toolkit's routine returns 999992.17.

## Score, information, and orthogonality

Writing \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\ for the derivative of
\\\log I_0\\, \$\$\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu),
\qquad \dfrac{\partial\ell}{\partial\kappa} = \cos(y-\mu) -
A(\kappa),\$\$ and every second derivative is closed form, with
\\A'(\kappa) = 1 - A(\kappa)/\kappa - A(\kappa)^2\\ from the Bessel
recurrences and no further evaluation. Since \\\mathbb{E}\[\sin(Y-\mu)\]
= 0\\ by symmetry, **the two parameters are orthogonal**: the expected
information is diagonal, and Fisher scoring updates the direction and
the concentration independently.

\\A(\kappa)\\ is the mean resultant length, so \\A'(\kappa)\\ is the
variance of \\\cos(Y-\mu)\\ and is positive, which keeps the information
positive definite.

The log-density is linear in \\\kappa\\ apart from the normalizing
constant, so at orders three and four every component naming one \\\mu\\
and two or more \\\kappa\\ is exactly zero.

## The moments a generic returns are not the circular ones

[`mean.distrib()`](https://statmodels7.github.io/distributions7/reference/mean.distrib.md),
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
are the ordinary moments of \\Y\\ as a number on \\\[-\pi, \pi)\\,
obtained numerically. They are **not** the circular quantities that
\\\mu\\ and \\\kappa\\ describe. \\\mu\\ is the mean *direction*, and
\\\mathbb{E}\[Y\] \ne \mu\\ whenever \\\mu \ne 0\\, because the interval
is cut at \\\pm\pi\\ rather than at \\\mu \pm \pi\\ and the density is
not symmetric about \\\mu\\ on it: at \\\mu = 1.2\\ and \\\kappa = 2\\
the ordinary mean is 1.079.

The circular mean is \\\mu\\ and the mean resultant length is \\\rho =
I_1(\kappa)/I_0(\kappa)\\, both closed form. Neither is returned by a
generic whose name means something else; compute them from a sample as
`atan2(mean(sin(z)), mean(cos(z)))` and
`sqrt(mean(cos(z))^2 + mean(sin(z))^2)`.

## The direction is carried on a bounded chart

The default link maps the free scale onto \\(-\pi, \pi)\\, which keeps
\\\mu\\ identified. The cost is that a fit cannot walk across the
boundary, so data concentrated near \\\pm\pi\\ are better rotated before
fitting than handed to the optimizer as they are. Leaving \\\mu\\
unbounded would make the likelihood periodic and every maximum one of
infinitely many.

## The distribution function

The density has no elementary antiderivative. The base class would
integrate it once per observation, which measured 4.5 seconds at ten
thousand points; this class supplies
[`distrib_cdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises1Distrib.md)
instead, a Fourier series integrated term by term whose term count is
\\8.5\sqrt{\kappa} + 10\\, measured. The **quantile** does come from
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md),
by root finding on that series.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\kappa \> 0\\ the concentration, \\I_m\\ the modified
Bessel function of the first kind of order \\m\\, and \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\ the mean resultant length. \\\eta\\ is a
parameter on the unconstrained scale of its link, with \\\theta =
g^{-1}(\eta)\\.

## References

Best, D. J. and Fisher, N. I. (1979). Efficient simulation of the von
Mises distribution. *Journal of the Royal Statistical Society, Series
C*, **28**(2), 152-157.

Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*, Chapter
3. Wiley, Chichester.

## See also

[`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
for the same law in the mean resultant length;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the analogous family on the line, which the von Mises approaches at
a large concentration;
[`numericals7::log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.html)
and
[`numericals7::bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.html)
for the Bessel machinery;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[VonMises1Distrib](https://statmodels7.github.io/distributions7/reference/VonMises1Distrib.md)
for the class.

## Examples

``` r
d <- vonmises1_distrib()
d
#> Distribution: Von Mises1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (mean direction)     | Link: bounded(lwr=-3.14159265358979, upr=3.14159265358979) | Domain: (-3.14159265358979, 3.14159265358979)
#>   kappa (concentration)      | Link: log        | Domain: (0, Inf)

# The density integrates to one over the circle.
th <- list(mu = 0.5, kappa = 2)
integrate(function(v) distrib_pdf(d, v, th), -pi, pi)$value
#> [1] 1

# The expected information is diagonal: direction and concentration are
# orthogonal, the sine having mean zero.
vapply(distrib_expected_hessian(d, 0, th), function(v) v[1], numeric(1))
#>       mu_mu    mu_kappa kappa_kappa 
#>  -1.3955493   0.0000000  -0.1642232 

# The ordinary mean is not the mean direction, the interval being cut at
# plus or minus pi.
c(ordinary = mean(d, list(mu = 1.2, kappa = 2)), direction = 1.2)
#>  ordinary direction 
#>  1.079473  1.200000 

# A sample recovers the circular quantities, which no generic returns.
set.seed(1)
z <- distrib_rng(d, 3e5, th)
c(circular_mean = atan2(mean(sin(z)), mean(cos(z))), mu = 0.5)
#> circular_mean            mu 
#>     0.4990003     0.5000000 
c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2),
  A = numericals7::bessel_i_ratio(2))
#> resultant         A 
#> 0.6982004 0.6977747 

# Fitting recovers both parameters.
set.seed(7)
coef(fit_distrib(d, distrib_rng(d, 2000, list(mu = 0.8, kappa = 3))))
#>        mu     kappa 
#> 0.7777677 3.0723246 
```
