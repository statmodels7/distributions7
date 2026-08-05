# von Mises Distribution Object

Creates a distribution object for the von Mises distribution, the
natural family for an angle, parametrised by a mean direction \\\mu\\
and a concentration \\\kappa\\.

## Usage

``` r
vonmises_distrib(
  link_mu = bounded_link(lwr = -pi, upr = pi),
  link_kappa = log_link()
)
```

## Arguments

- link_mu:

  A link function object for \\\mu\\. Defaults to
  [`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.html)
  on \\(-\pi, \pi)\\.

- link_kappa:

  A link function object for \\\kappa\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `VonMisesDistrib`.

## Details

The observation is an angle and the support is a circle, written here as
\\\[-\pi, \pi)\\. This is the first family in the package whose support
has that topology: the two ends of the interval are the same point, so a
density need not vanish at either.

**Density:** \$\$f(y) = \dfrac{e^{\kappa\cos(y-\mu)}}{2\pi
I_0(\kappa)}\$\$ The normalising constant is a modified Bessel function,
and it is evaluated exponentially scaled with the exponent added back,
so a concentration past \\\kappa = 700\\ does not overflow.

**Score, observed and expected Hessian.** Writing \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\ for the derivative of \\\log I_0\\,
\$\$\dfrac{\partial\ell}{\partial\mu} = \kappa\sin(y-\mu), \qquad
\dfrac{\partial\ell}{\partial\kappa} = \cos(y-\mu) - A(\kappa),\$\$ and
every second derivative is closed form, with \\A'(\kappa) = 1 -
A(\kappa)/\kappa - A(\kappa)^2\\ obtained from the Bessel recurrences
rather than from a further evaluation. Since \\\mathbb{E}\[\sin(Y-\mu)\]
= 0\\ by symmetry, **the two parameters are orthogonal**: the expected
information is diagonal, and Fisher scoring updates the direction and
the concentration independently.

\\A(\kappa)\\ is the mean resultant length, so \\A'(\kappa)\\ is the
variance of \\\cos(Y-\mu)\\ and is positive, which is what makes the
information positive definite.

**The mean direction is carried on a bounded chart**, the default link
mapping the free scale onto \\(-\pi, \pi)\\. That keeps the parameter
identified, at the cost that a fit cannot walk across the boundary: data
concentrated near \\\pm\pi\\ are better rotated before fitting than
handed to the optimiser as they are. Leaving \\\mu\\ unbounded instead
would make the likelihood periodic and every maximum one of infinitely
many.

**Parameter domains:**

- \\\mu \in (-\pi, \pi)\\

- \\\kappa \in (0, +\infty)\\

The distribution function has no elementary form and comes from the base
class by quadrature over the bounded support, with the quantile by root
finding on it.

## References

Best, D. J. and Fisher, N. I. (1979). Efficient simulation of the von
Mises distribution. *Applied Statistics* 28, 152-157.

Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*. Wiley.

## See also

[`gaussian_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian_distrib.md)
for the analogous family on the line

## Examples

``` r
d <- vonmises_distrib()
d@params
#> [1] "mu"    "kappa"

theta <- list(mu = 0.5, kappa = 2)
distrib_pdf(d, c(-1, 0, 0.5, 2), theta)
#> [1] 0.08042773 0.40385253 0.51588541 0.08042773

# the expected information is diagonal: direction and concentration are
# orthogonal, the sine having mean zero
distrib_expected_hessian(d, 0, theta)
#> $mu_mu
#> [1] -1.395549
#> 
#> $mu_kappa
#> [1] 0
#> 
#> $kappa_kappa
#> [1] -0.1642232
#> 
```
