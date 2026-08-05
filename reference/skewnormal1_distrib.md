# Skew Normal Distribution Object

Creates a distribution object for Azzalini's skew normal distribution,
with location \\\mu\\, scale \\\sigma\\ and shape \\\alpha\\. The
gaussian is the special case \\\alpha = 0\\.

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

  A link function object for the location \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link function object for the scale \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

- link_alpha:

  A link function object for the shape \\\alpha\\, which is
  unconstrained. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

## Value

An S7 object of class
[`SkewNormal1Distrib`](https://statmodels7.github.io/distributions7/reference/SkewNormal1Distrib.md)
(inheriting from `continuous_distrib`).

## Details

**Probability density function**, with \\z = (y-\mu)/\sigma\\: \$\$f(y;
\mu, \sigma, \alpha) = \dfrac{2}{\sigma}\\\phi(z)\\\Phi(\alpha z)\$\$
The factor \\2\Phi(\alpha z)\\ tilts the gaussian: it is above one where
\\\alpha z \> 0\\ and below one where \\\alpha z \< 0\\, so positive
\\\alpha\\ skews the density to the right. At \\\alpha = 0\\ it is
identically one and the family reduces to the gaussian.

**Cumulative distribution function:** \$\$F(q; \mu, \sigma, \alpha) =
\Phi(z) - 2\\T(z, \alpha)\$\$ with \\T\\ Owen's T function; see
[`owen_t`](https://statmodels7.github.io/numericals7/reference/owen_t.html).
The quantile function has no closed form and comes from the base class
by root finding.

**Score and observed Hessian** are closed form, written in the inverse
Mills ratio \\R(t) = \phi(t)/\Phi(t)\\ at \\t = \alpha z\\ and its
derivative \\R' = -R(t+R)\\; see
[`distrib_gradient.SkewNormal1Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md)
and
[`distrib_hessian.SkewNormal1Distrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal1Distrib.md).

**Expected information.** There is no elementary closed form, so none is
registered and
[`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
approximates it by the strategy named in `approx`, the default being the
score variance.

**Singularity at the gaussian.** At \\\alpha = 0\\ the expected
information of this parametrisation is **singular**: the derivative in
\\\alpha\\ becomes collinear with the derivative in \\\mu\\, so the two
cannot be separated there. This is a property of the family and not of
the implementation, and it is why the profile log-likelihood in
\\\alpha\\ is flat at the origin. A fit whose true shape is near zero
will report a large standard error for \\\alpha\\; the centred
parametrisation of Azzalini and Capitanio removes the singularity and is
a different object, not a reparametrisation this class performs.

**Moments.** With \\\delta = \alpha/\sqrt{1+\alpha^2}\\ and \\b =
\sqrt{2/\pi}\\, the mean is \\\mu + \sigma b \delta\\ and the variance
\\\sigma^2(1 - b^2\delta^2)\\. The skewness is bounded: it lies in
\\(-0.9953, 0.9953)\\ whatever \\\alpha\\ is, which is the limitation of
the family and the reason the skew \\t\\ exists.

**Higher orders.** The observed third and fourth derivatives are closed
form, every derivative of \\\log\Phi(t)\\ being a polynomial in \\t\\
and the inverse Mills ratio through \\R' = -R(t+R)\\; their expected
values share the obstruction of the expected information and are
approximated numerically.

**Parameter Domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

- \\\alpha \in (-\infty, +\infty)\\

## References

Azzalini, A. (1985). A class of distributions which includes the normal
ones. *Scandinavian Journal of Statistics* 12, 171-178.

Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
Families*. Cambridge University Press.

## Examples

``` r
d <- skewnormal1_distrib()
d@params
#> [1] "mu"    "sigma" "alpha"

theta <- list(mu = 0, sigma = 1, alpha = 3)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.0006532716 0.3989422804 0.4832881774
distrib_gradient(d, c(-1, 0, 1), theta)
#> $mu
#> [1] -10.8492960  -2.3936537   0.9866865
#> 
#> $sigma
#> [1]  9.84929596 -1.00000000 -0.01331352
#> 
#> $alpha
#> [1] -3.283098655  0.000000000  0.004437839
#> 

# shape zero is the gaussian
max(abs(distrib_pdf(d, c(-1, 0, 1), list(mu = 0, sigma = 1, alpha = 0)) -
        stats::dnorm(c(-1, 0, 1))))
#> [1] 5.551115e-17

# the skewness the family can reach is bounded
c(alpha_3 = skewness(d, theta),
  alpha_50 = skewness(d, list(mu = 0, sigma = 1, alpha = 50)))
#>   alpha_3  alpha_50 
#> 0.6670236 0.9936306 
```
