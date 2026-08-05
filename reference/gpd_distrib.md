# Generalized Pareto Distribution Object

Creates a distribution object for the generalized Pareto distribution,
parametrized by a scale \\\sigma\\ and a shape \\\xi\\.

## Usage

``` r
gpd_distrib(link_sigma = log_link(), link_xi = identity_link())
```

## Arguments

- link_sigma:

  A link function object for \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

- link_xi:

  A link function object for \\\xi\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the shape being free to take either sign.

## Value

An S7 object of class `GPDDistrib`.

## Details

The family of exceedances over a high threshold, and the natural
companion of
[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
in an analysis of extremes.

**Density and distribution function:** \$\$f(y) =
\dfrac{1}{\sigma}\left(1 + \dfrac{\xi y}{\sigma}\right)^{-1/\xi-1},
\qquad F(q) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi}\$\$ At
\\\xi = 0\\ both reduce to the exponential, which the implementation
reaches through a series rather than by a special case, so the parameter
may pass through zero during a fit.

**It is not parametrized by its mean**, unlike most families here. The
mean is \\\sigma/(1-\xi)\\ and exists only for \\\xi \< 1\\, so a mean
parametrization would leave the family undescribable exactly where it is
most used — the heavy-tailed regime. This is the argument that keeps
[`mv_sigma`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
and
[`variance`](https://statmodels7.github.io/distributions7/reference/variance.md)
apart for the multivariate \\t\\: a parametrization must not depend on a
moment that need not exist.

**The support depends on the parameters when \\\xi \< 0\\**, being
\\\[0, -\sigma/\xi\]\\, and this is the first family here of which that
is true. What it costs is the automatic license to differentiate under
the integral sign, on which the Bartlett identities rest. Two things
survive and one does not:

- the derivatives returned are correct as derivatives of the log-density
  at every admissible point, whatever the sign of \\\xi\\;

- the expected information exists and is the closed form above for \\\xi
  \> -1/2\\. The condition is exactly that the integrand be integrable:
  near the upper endpoint the second derivative grows like
  \\(1-u)^{-2\|\xi\|}\\ in the probability scale, which is integrable if
  and only if \\\|\xi\| \< 1/2\\;

- below \\\xi = -1/2\\ the information does not exist,
  [`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
  returns `NA`, and the classical asymptotics of the maximum likelihood
  estimator do not hold (Smith, 1985).

The `bounds` of the object are `c(0, Inf)` because they are fixed at
construction while the true endpoint moves with the parameters; the
density is zero beyond it, so nothing computes a wrong number, but a
caller reading `bounds` learns less than usual.

**Moments:** mean \\\sigma/(1-\xi)\\ for \\\xi \< 1\\, variance
\\\sigma^2/\\(1-\xi)^2(1-2\xi)\\\\ for \\\xi \< 1/2\\.

**Parameter domains:**

- \\\sigma \in (0, +\infty)\\

- \\\xi \in (-\infty, +\infty)\\

## References

Smith, R. L. (1985). Maximum likelihood estimation in a class of
nonregular cases. *Biometrika* 72, 67-90.

Davison, A. C. and Smith, R. L. (1990). Models for exceedances over high
thresholds. *Journal of the Royal Statistical Society B* 52, 393-442.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md),
[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md),
[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)

## Examples

``` r
d <- gpd_distrib()
d@params
#> [1] "sigma" "xi"   

theta <- list(sigma = 1.5, xi = 0.3)
distrib_pdf(d, c(0.2, 1, 4), theta)
#> [1] 0.5624677 0.3025450 0.0522069
distrib_gradient(d, c(0.2, 1, 4), theta)
#> $sigma
#> [1] -0.5555556 -0.1851852  0.6172840
#> 
#> $xi
#> [1] -0.1197699 -0.3816123  0.1112099
#> 

# at xi = 0 it is the exponential, reached by a series
max(abs(distrib_pdf(d, c(0.2, 1, 4), list(sigma = 1.5, xi = 0)) -
        dexp(c(0.2, 1, 4), 1 / 1.5)))
#> [1] 1.387779e-17

# the information exists only above -1/2
distrib_expected_hessian(d, 0, list(sigma = 1.5, xi = -0.7))
#> $sigma_sigma
#> [1] NA
#> 
#> $sigma_xi
#> [1] NA
#> 
#> $xi_xi
#> [1] NA
#> 
```
