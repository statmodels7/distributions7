# Skew Normal Distribution Object, Centered Parametrization

Builds a skew normal distribution object parametrized by its mean
\\\mu\\, its standard deviation \\\sigma\\ and its skewness
\\\gamma_1\\. It is the same family as
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
named by three moments instead of by a location, a scale and a shape.

The parametrization is Azzalini's. Its expected information stays
non-singular at \\\gamma_1 = 0\\, where the direct parametrization's
loses a rank. The map between the two runs through a cube root, so the
derivative of the shape in \\\gamma_1\\ is unbounded at zero skewness
and the log-density has no derivative of order two or more in
\\\gamma_1\\ there.

## Usage

``` r
skewnormal2_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_gamma1 = bounded_link(lwr = -sn_max_skew(), upr = sn_max_skew())
)
```

## Arguments

- link_mu:

  A `linkfunctions7` link object for the mean, which is unconstrained.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link object for the standard deviation, which must be strictly
  positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_gamma1:

  A link object for the skewness, which must lie strictly inside
  \\(-0.9952717, 0.9952717)\\. Defaults to
  [`linkfunctions7::bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.html)
  over exactly that interval, so any real linear predictor maps to a
  skewness the family can reach.

## Value

An S7 object of class
[SkewNormal2Distrib](https://statmodels7.github.io/distributions7/reference/SkewNormal2Distrib.md),
inheriting from `continuous_distrib`. Its `params` are
`c("mu", "sigma", "gamma1")`, its `bounds` `c(-Inf, Inf)`, and its
`link_params` the three links given here.

## Relation to reparametrize()

The map passes through \\c =
\mathrm{sign}(\gamma_1)(2\|\gamma_1\|/(4-\pi))^{1/3}\\, and two things
follow from it. It carries a sign, so
[`sn2_theta()`](https://statmodels7.github.io/distributions7/reference/sn2_theta.md)
takes the sign from the plain value of `gamma1` and passes it to
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
as an argument, and the body of the map contains no
[`abs()`](https://rdrr.io/r/base/MathFun.html). The derivatives of the
log-density in the centered parameters are written out per order in
compiled kernels (see
[`distrib_gradient.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md));
near \\\gamma_1 = 0\\, where their closed forms cancel, they come from
the series of the log-density in \\\gamma_1^{1/3}\\, and far below the
step of the skewing function from a composition that does not cancel.

## Behavior near zero skewness

\\\partial\alpha/\partial\gamma_1\\ grows without bound as \\\gamma_1
\to 0\\, like \\\gamma_1^{-2/3}\\. The score does not follow it: the
divergent contributions cancel, so the information in \\\gamma_1\\ tends
to \\1/6\\ and the whole matrix stays positive definite at symmetry.

The **observed** curvature does diverge, at the rate the cube root sets:
\\\gamma_1^{-2/3}\\. The derivatives of order two and more are therefore
rejected at \\\gamma_1 = 0\\ exactly, with a message naming the cause;
the score and the expected information are returned there. The density,
the distribution function, the quantile function, the generator and both
response derivatives are fine there and equal the Gaussian's.

For \\\|\gamma_1\|\\ below
[`sn2_ge()`](https://statmodels7.github.io/distributions7/reference/sn2_ge.md)
the expected information is computed from series in \\\gamma_1^{1/3}\\,
which contain no cancellation. A problem that is symmetric by
construction is better posed in
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
since the observed derivatives of order two and more in the centered
parametrization are rejected at \\\gamma_1 = 0\\.

## The bound on the skewness

A skew normal cannot reach \\\|\gamma_1\| \> 0.9952717\\ whatever its
shape, so `gamma1` is bounded there and carries a bounded link by
default. A model that needs a larger skewness can use
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

## Parameter domains

- \\\mu \in (-\infty, \infty)\\

- \\\sigma \in (0, \infty)\\

- \\\gamma_1 \in (-0.9952717, 0.9952717)\\

## The distribution

\$\$f(y) =
\frac{2}{\omega}\\\phi\\\left(\frac{y-\xi}{\omega}\right)\Phi\\\left(\alpha\\\frac{y-\xi}{\omega}\right),
\qquad (\xi, \omega, \alpha) = \mathrm{DP}(\mu, \sigma, \gamma_1)\$\$ on
\\y \in \mathbb{R}\\, with \\\mathrm{DP}\\ the map of
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md).

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2},
\qquad \gamma_1(Y) = \gamma_1.\$\$

## References

Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
Families*. Cambridge University Press. The centered parametrization is
section 3.1.4.

## See also

[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the direct parametrization,
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
for the map between them,
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
for a family that reaches a larger skewness, and
[SkewNormal2Distrib](https://statmodels7.github.io/distributions7/reference/SkewNormal2Distrib.md)
for the class and its method list.

## Examples

``` r
d <- skewnormal2_distrib()
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# All three parameters are moments, which is what "centered" names.
c(mean = mean(d, th), sd = sqrt(variance(d, th)),
  skewness = skewness(d, th))
#>     mean       sd skewness 
#>      0.0      1.0      0.5 

# The fourth moment is not free: it follows from the skewness.
vapply(c(0.001, 0.5, 0.9),
       function(g) kurtosis(d, list(mu = 0, sigma = 1, gamma1 = g)), 0)
#> [1] 8.746873e-05 3.471199e-01 7.600512e-01

# The information stays invertible at symmetry, where the direct
# parametrization's does not.
smallest_eigenvalue <- function(dd, p) {
  e <- distrib_expected_hessian(dd, 0, p)
  M <- matrix(c(e[[1]], e[[4]], e[[5]],
                e[[4]], e[[2]], e[[6]],
                e[[5]], e[[6]], e[[3]]), 3, 3)   # hess_names() order
  min(eigen(-M, only.values = TRUE)$values)
}
c(centered = smallest_eigenvalue(d, list(mu = 0, sigma = 1, gamma1 = 1e-6)),
  direct = smallest_eigenvalue(skewnormal1_distrib(),
                               list(mu = 0, sigma = 1, alpha = 0)))
#>      centered        direct 
#>  1.666723e-01 -2.550320e-27 

# A fit recovers all three moments.
set.seed(11)
x <- distrib_rng(d, 4000, list(mu = 3, sigma = 2, gamma1 = 0.6))
coef(fit_distrib(d, x))
#>        mu     sigma    gamma1 
#> 2.9797006 1.9525233 0.5993477 
```
