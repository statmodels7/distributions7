# von Mises Distribution, Mean Direction and Mean Resultant Length

Builds the distribution object for the von Mises family parametrized by
its mean direction \\\mu\\ and its **mean resultant length** \\\rho \in
(0, 1)\\. The returned object carries closed-form derivatives of the
log-density to fourth order, a closed-form expected information, and the
same Bessel-series distribution function the concentration
parametrization uses.

The concentration \\\kappa\\ of
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
is unbounded and hard to read. The resultant length is bounded, is the
quantity circular statistics reports, and is one minus the circular
variance; a sample reads it back directly as
`sqrt(mean(cos(z))^2 + mean(sin(z))^2)`.

## Usage

``` r
vonmises2_distrib(
  link_mu = bounded_link(lwr = -pi, upr = pi),
  link_rho = logit_link()
)
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean direction \\\mu\\.
  Defaults to `linkfunctions7::bounded_link(lwr = -pi, upr = pi)`. The
  chart keeps \\\mu\\ identified at the cost that a fit cannot walk
  across \\\pm\pi\\; see
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- link_rho:

  A `link` object from `linkfunctions7` for the resultant length
  \\\rho\\. Defaults to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  the natural link onto \\(0, 1)\\.

## Value

An S7 object of class `VonMises2Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"von mises2"`, `dimension`
`"univariate"`, `bounds` `c(-pi, pi)`, `params` `c("mu", "rho")`,
`n_params` `2`, `params_bounds` the domains \\(-\pi, \pi)\\ and \\(0,
1)\\, and `link_params` the two links given here.

## The parametrization, and why it is a family of its own

The density on \\y \in \[-\pi, \pi)\\ is \$\$f(y; \mu, \rho) =
\frac{e^{\kappa\cos(y-\mu)}}{2\pi I\_{0}(\kappa)}, \qquad \kappa =
A^{-1}(\rho),\$\$ with \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\, and
\\\mathbb{E}\[\cos(Y-\mu)\] = \rho\\. The map \\\rho = A(\kappa)\\ is a
strictly increasing bijection from \\(0, \infty)\\ onto \\(0, 1)\\.

**Its inverse has no closed form**, which is why this is a family of its
own and not a
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
of the other. \\\kappa\\ is obtained by root finding on \\\log\kappa\\,
and its four derivatives come from the inverse function rule applied to
\\A'\\ through \\A^{(4)}\\, which the Bessel recurrences give from the
same two evaluations \\A\\ already needs.

The map is steep near \\\rho = 1\\: measured,
\\\mathrm{d}\kappa/\mathrm{d}\rho\\ is 2.00 at \\\rho = 0.01\\, 3.14 at
0.5, 199 at 0.95 and \\5.0\times10^{5}\\ at 0.999. A fit whose data are
nearly all in one direction is therefore far better conditioned in
\\\kappa\\ than in \\\rho\\.

## Derivatives and information

The map touches the **second parameter only**, so every chain rule is
the one-variable one and the derivatives are exact at every order. Every
component carrying at least one \\\mu\\ collapses to a single term, the
concentration parametrization's \\\mu\\-derivatives being linear in
\\\kappa\\.

The expected information is closed form and the two parameters are
orthogonal, as in the concentration parametrization, with
\$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\kappa\rho, \qquad
\mathbb{E}\[\ell^{(\mu\rho)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\rho\rho)}\] = -\dfrac{1}{A'(\kappa)},\$\$ the last
being the reciprocal of the information in \\\kappa\\, since the
Jacobian of the map is that reciprocal.

## The moments are not the parameters

[`mean.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.VonMises2Distrib.md)
returns the ordinary expectation of \\Y\\ on \\\[-\pi, \pi)\\, which
differs from \\\mu\\ whenever \\\mu \ne 0\\;
[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
discusses the same distinction. What this parametrization buys is that
the **second** parameter is a quantity a sample reads back directly,
which the concentration is not.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\rho \in (0,1)\\ the mean resultant length, \\\kappa\\ the
concentration, \\I_m\\ the modified Bessel function of the first kind of
order \\m\\, and \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\. \\\eta\\ is a
parameter on the unconstrained scale of its link, with \\\theta =
g^{-1}(\eta)\\.

## References

Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*, Chapter
3. Wiley, Chichester.

## See also

[`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
for the same law in the concentration, which is better conditioned at a
nearly deterministic direction;
[`numericals7::bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.html)
for the map and its derivatives;
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for another family on \\(0, 1)\\;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[VonMises2Distrib](https://statmodels7.github.io/distributions7/reference/VonMises2Distrib.md)
for the class.

## Examples

``` r
d <- vonmises2_distrib()
d
#> Distribution: Von Mises2
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu  (mean direction)     | Link: bounded(lwr=-3.14159265358979, upr=3.14159265358979) | Domain: (-3.14159265358979, 3.14159265358979)
#>   rho (mean resultant length) | Link: logit      | Domain: (0, 1)

# The resultant length is bounded, which is what makes it readable.
d@params_bounds$rho
#> [1] 0 1

# The same law as the concentration parametrization at the implied kappa.
th <- list(mu = 0.5, rho = 0.7)
k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
all.equal(distrib_pdf(d, c(-1, 0, 1), th),
          distrib_pdf(vonmises1_distrib(), c(-1, 0, 1),
                      list(mu = 0.5, kappa = k)))
#> [1] TRUE

# A sample reads the second parameter back directly, which is the point.
set.seed(1)
z <- distrib_rng(d, 3e5, th)
c(resultant = sqrt(mean(cos(z))^2 + mean(sin(z))^2), rho = 0.7)
#> resultant       rho 
#> 0.7003239 0.7000000 

# The map is steep near one, so a nearly deterministic direction is better
# conditioned in kappa than in rho.
vapply(c(0.01, 0.5, 0.95, 0.999),
       function(r) numericals7::bessel_i_ratio_inverse(r)$d1, numeric(1))
#> [1] 2.000300e+00 3.137622e+00 1.994903e+02 4.999996e+05

# Fitting recovers both parameters.
set.seed(3)
coef(fit_distrib(d, distrib_rng(d, 2000, list(mu = 0.8, rho = 0.6))))
#>        mu       rho 
#> 0.7674996 0.5837115 
```
