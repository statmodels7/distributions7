# Skew Normal Distribution Object, Centered Parametrization

Builds a skew normal distribution object parametrized by its mean
\\\mu\\, its standard deviation \\\sigma\\ and its skewness
\\\gamma_1\\. It is the same family as
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
named by three moments instead of by a location, a scale and a shape.

The parametrization is Azzalini's, and the property it exists for is
that its expected information stays non-singular at \\\gamma_1 = 0\\,
where the direct parametrization's loses a rank. The price is that the
map between the two runs through a cube root, so no parameter derivative
exists at exactly zero skewness.

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

## This is a family, not a reparametrize()

The map passes through \\c =
\mathrm{sign}(\gamma_1)(2\|\gamma_1\|/(4-\pi))^{1/3}\\, and two things
follow from it. It carries a sign, so
[`sn2_theta()`](https://statmodels7.github.io/distributions7/reference/sn2_theta.md)
reads that off the plain value and hands it to
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
as an argument, leaving a body with no
[`abs()`](https://rdrr.io/r/base/MathFun.html) in it to differentiate.
And its derivatives are written out by hand in
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md),
as a keyed table of the map's partials, which
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
consumes. The toolkit assembles higher-order derivatives from
written-out tables throughout.

## What the map costs, and what it buys

\\\partial\alpha/\partial\gamma_1\\ grows without bound as \\\gamma_1
\to 0\\: measured, 3.9 at \\\gamma_1 = 0.5\\, 12.8 at 0.01 and 258 at
\\10^{-4}\\. The score does not follow it: the divergent contributions
cancel, so the information in \\\gamma_1\\ tends to \\1/6\\ and the
whole matrix stays positive definite at symmetry.

The **observed** curvature does diverge, at the rate the cube root sets:
\\\gamma_1^{-2/3}\\, measured at 4.642 per decade against \\10^{2/3} =
4.6416\\. The parameter derivatives are therefore rejected at \\\gamma_1
= 0\\ exactly, with a message naming the cause. The density, the
distribution function, the quantile function, the generator and both
response derivatives are fine there and equal the Gaussian's.

The cancellation is between terms of size \\\gamma_1^{-2/3}\\, so it
eventually runs out of digits. Measured, the expected information's
\\\gamma_1\\ component holds seven figures at \\\gamma_1 = 10^{-8}\\,
loses three by \\10^{-10}\\ and is negative at \\10^{-12}\\. Those are
values no fit visits, and a genuinely symmetric problem is better posed
in
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## The bound on the skewness

A skew normal cannot reach \\\|\gamma_1\| \> 0.9952717\\ whatever its
shape, so `gamma1` is bounded there and carries a bounded link by
default. That ceiling is why
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
exists.

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
#>  1.666723e-01 -5.738220e-27 

# A fit recovers all three moments.
set.seed(11)
x <- distrib_rng(d, 4000, list(mu = 3, sigma = 2, gamma1 = 0.6))
coef(fit_distrib(d, x))
#>        mu     sigma    gamma1 
#> 2.9797006 1.9525233 0.5993477 
```
