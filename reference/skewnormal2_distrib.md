# Skew Normal Distribution in Its Centered Parametrization

Creates a skew normal distribution object parametrized by its mean, its
standard deviation and its skewness.

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

  Link function for the mean. Defaults to the identity.

- link_sigma:

  Link function for the standard deviation. Defaults to the log.

- link_gamma1:

  Link function for the skewness. Defaults to a link bounded to the
  reachable interval.

## Value

An S7 object of class
[`SkewNormal2Distrib`](https://statmodels7.github.io/distributions7/reference/SkewNormal2Distrib.md).

## Details

The direct parametrization of
[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
carries a location, a scale and a shape, none of which is a moment. Here
all three parameters are moments, and the family is Azzalini's centered
parametrization.

**Why this is not a reparametrize().** The map passes through \\c =
\mathrm{sign}(\gamma_1)(2\|\gamma_1\|/(4-\pi))^{1/3}\\, and two things
follow. It carries a sign, which a jet cannot take of itself, so the
sign is read off the plain value before any jet is seeded. And
\\\partial\alpha/\partial\gamma_1\\ grows without bound as \\\gamma_1
\to 0\\: measured, 3.9 at \\\gamma_1 = 0.5\\ and 258 at \\10^{-4}\\. The
value of the parametrization is that the score does **not** follow it,
the divergent contributions canceling, so the variance of the score in
\\\gamma_1\\ is 0.158 at \\\gamma_1 = 0.05\\ and 0.158 again at
\\0.01\\.

**What that costs in arithmetic.** The cancellation is between terms of
size proportional to the Jacobian, so the significant digits lost grow
like the logarithm of it: negligible over the range a fit visits, and
severe only within a few multiples of \\10^{-8}\\ of exact symmetry. The
**expected information**, unlike the direct parametrization's, is
non-singular at zero skewness, which is the property the parametrization
exists for.

**The bound on the skewness.** A skew normal cannot reach \\\|\gamma_1\|
\> 0.9952717\\, whatever \\\alpha\\ is, so the parameter is bounded
there. That ceiling is the reason the skew \\t\\ exists.

## The distribution

\$\$f(y) =
\frac{2}{\omega}\\\phi\\\left(\frac{y-\xi}{\omega}\right)\Phi\\\left(\alpha\\\frac{y-\xi}{\omega}\right),
\qquad (\xi, \omega, \alpha) = \mathrm{DP}(\mu, \sigma, \gamma_1)\$\$ on
\\y \in \mathbb{R}\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2}\$\$

## References

Azzalini, A. and Capitanio, A. (2014). *The Skew-Normal and Related
Families*. Cambridge University Press. The centred parametrization is
section 3.1.4.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)

## Examples

``` r
d <- skewnormal2_distrib()
theta <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# all three parameters are moments, which is what centered means
c(mean = mean(d, theta), sd = sqrt(variance(d, theta)),
  skewness = skewness(d, theta))
#>     mean       sd skewness 
#>      0.0      1.0      0.5 
```
