# Pseudo-Huber Distribution, Location, Scale and Shape

Builds the distribution object for the pseudo-Huber family, whose
log-density is the negative pseudo-Huber loss \\-\sqrt{\nu +
\\(y-\mu)/\sigma\\^2}\\. It is the likelihood counterpart of that loss:
the score is bounded, so a gross outlier contributes a limited amount to
the estimating equation, and the shape \\\nu \> 0\\ interpolates between
a Laplace at \\\nu \to 0\\ and a Gaussian at \\\nu \to \infty\\. The
family is the symmetric hyperbolic distribution, a special case of the
generalized hyperbolic.

The returned object carries closed-form derivatives of the log-density
to fourth order in the parameters and closed first and second
derivatives in the response. The distribution function, the quantile and
the expected information are numerical; that is the price of the family.

## Usage

``` r
pseudohuber_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_nu = log_link()
)
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the location \\\mu\\.
  Defaults to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the location ranging over the whole line already.

- link_sigma:

  A `link` object from `linkfunctions7` for the scale \\\sigma\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

- link_nu:

  A `link` object from `linkfunctions7` for the shape \\\nu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `PseudoHuberDistrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"pseudo huber"`, `dimension`
`"univariate"`, `bounds` `c(-Inf, Inf)`, `params`
`c("mu", "sigma", "nu")`, `n_params` `3`, `params_bounds` the domains
\\(-\infty, \infty)\\, \\(0, \infty)\\ and \\(0, \infty)\\, and
`link_params` the three links given here.

## The parametrization

The density on \\y \in (-\infty, \infty)\\ is \$\$f(y; \mu, \sigma, \nu)
= \dfrac{1}{2 \sigma \sqrt{\nu}\\ K_1(\sqrt{\nu})}
\exp\left(-\sqrt{\nu + \left(\dfrac{y-\mu}{\sigma}\right)^2}\right),\$\$
with \\K_1\\ the modified Bessel function of the second kind. The
exponent is quadratic in the residual near the location, where
\\\sqrt{\nu + z^2} \approx \sqrt{\nu} + z^2/(2\sqrt{\nu})\\, and linear
far from it, where it is \\\|z\|\\. That is the pseudo-Huber loss, and
this family is its exponential.

The Bessel terms are degree-homogeneous, so the exponentially scaled
`besselK(x, nu, expon.scaled = TRUE)` is exact and avoids the overflow
the unscaled form meets; the constant stays finite to \\\nu = 2000\\.

## The two limits, with their rates

At \\\nu \to 0\\ the exponent is \\\|y-\mu\|/\sigma\\ and the family is
the **Laplace** with scale \\\sigma\\: measured at \\\mu = 0\\, \\\sigma
= 1\\, the largest gap over \\y \in \\0.5, 1, 2, 4\\\\ is `4.9e-05` at
\\\nu = 10^{-4}\\ and `1.9e-12` at \\10^{-12}\\, so the approach is
\\O(\nu)\\.

At \\\nu \to \infty\\ the exponent is \\\sqrt{\nu} +
(y-\mu)^2/(2\sqrt{\nu}\sigma^2)\\ and the family is the **Gaussian**
with standard deviation \\\sigma\nu^{1/4}\\: the largest gap is
`1.5e-04` at \\\nu = 10^4\\ and `1.5e-07` at \\10^8\\. Note that the
spread grows with \\\nu\\, so \\\sigma\\ alone is not the standard
deviation.

## Moments

Every moment is a ratio of Bessel functions at \\\sqrt{\nu}\\:
\$\$\mathbb{E}(Y) = \mu, \qquad \operatorname{Var}(Y) = \sigma^2
\sqrt{\nu}\\ \dfrac{K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})},\$\$ the skewness
is 0 by symmetry, and the excess kurtosis is
\\3K_3(\sqrt{\nu})K_1(\sqrt{\nu})/K_2(\sqrt{\nu})^2 - 3\\. The kurtosis
falls as \\\nu\\ grows, from 2.93 at \\\nu = 0.01\\ through 1.86 at 1 to
0.030 at \\10^4\\, tending to the Gaussian's 0.

## Derivatives

Write \\r = y - \mu\\ and \\D = \sqrt{\nu + (r/\sigma)^2}\\. The score
is \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2 D},
\qquad \dfrac{\partial \ell}{\partial \sigma} =
\dfrac{1}{\sigma}\left(\dfrac{r^2}{\sigma^2 D} - 1\right),\$\$
\$\$\dfrac{\partial \ell}{\partial \nu} = -\dfrac{1}{2}\left\[
\dfrac{1}{\nu} + \dfrac{1}{D} + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\\
K_1(\sqrt{\nu})}\right\].\$\$ The location component is **bounded by**
\\1/\sigma\\, which is the family's whole point; the curvature in
\\\mu\\, unlike a Student t's, is negative at every observation, so the
log-density stays concave in the location. Orders three and four are
closed form too, in compiled kernels, as are the derivatives in the
response and the mixed derivative
[`distrib_cross_y.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.PseudoHuberDistrib.md),
whose \\\nu\\ component is \\r/(2\sigma^2 D^3)\\.

## What is numerical, and what it costs

The distribution function has no elementary form and is a quadrature,
batched over the quantiles and reflected about \\\mu\\ so that only the
lower tail is integrated. The quantile inverts it by root-finding, and
the generator inverts it at uniform variates, so a sample costs one
root-find per draw.

The **expected information has no closed form either**. Its four
non-zero components come from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
and the two containing \\\mu\\ an odd number of times are replaced by
exact zeros, the law being symmetric. That makes the location orthogonal
to the scale and the shape.
[`expected_hessian_exact.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.PseudoHuberDistrib.md)
declares the approximation, so a caller who branches on the distinction
branches correctly; at 100 observations it costs about 11 seconds
against a median of 0.183 milliseconds for a family that writes its
information out.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Given the cost of the
expected information, `method = optimizers7::newton()` on the observed
Hessian is the cheaper route here, and it is closed form.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location,
\\\sigma \> 0\\ the scale, \\\nu \> 0\\ the shape, \\r = y - \mu\\, \\D
= \sqrt{\nu + (r/\sigma)^2}\\ and \\K_m\\ the modified Bessel function
of the second kind of order \\m\\, `besselK(x, m)` in R. \\\eta\\ is a
parameter on the unconstrained scale of its link, with \\\theta =
g^{-1}(\eta)\\.

## References

Barndorff-Nielsen, O. (1978). Hyperbolic distributions and distributions
on hyperbolae. *Scandinavian Journal of Statistics*, **5**(3), 151-157.

Charbonnier, P., Blanc-Feraud, L., Aubert, G. and Barlaud, M. (1997).
Deterministic edge-preserving regularization in computed imaging. *IEEE
Transactions on Image Processing*, **6**(2), 298-311.

## See also

[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the two limits;
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
for the other robust three-parameter family here, whose score redescends
to zero where this one flattens to a bound;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[PseudoHuberDistrib](https://statmodels7.github.io/distributions7/reference/PseudoHuberDistrib.md)
for the class.

## Examples

``` r
d <- pseudohuber_distrib()
d
#> Distribution: Pseudo Huber
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (location)           | Link: identity   | Domain: (-Inf, Inf)
#>   sigma (scale)              | Link: log        | Domain: (0, Inf)
#>   nu    (shape)              | Link: log        | Domain: (0, Inf)

# The density integrates to one.
th <- list(mu = 0.4, sigma = 1.2, nu = 2)
integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value
#> [1] 1

# Every moment is a Bessel ratio; the skewness is zero by symmetry.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>     mean      var     skew     kurt 
#> 0.400000 4.429997 0.000000 1.534651 

# A small shape is a Laplace of scale sigma; a large one a Gaussian of
# standard deviation sigma * nu^(1/4).
yy <- c(0.5, 1, 2, 4)
c(laplace = max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e-8)) -
                    0.5 * exp(-abs(yy)))),
  gaussian = max(abs(distrib_pdf(d, yy, list(mu = 0, sigma = 1, nu = 1e8)) -
                     dnorm(yy, 0, 100))))
#>      laplace     gaussian 
#> 1.186719e-08 1.495912e-07 

# The location score is bounded, where a Gaussian's grows without bound.
rr <- c(1, 4, 16, 64)
rbind(residual = rr,
      pseudohuber = distrib_gradient(d, 0.4 + rr, th)$mu,
      gaussian = rr / 1.2^2)
#>                  [,1]      [,2]      [,3]       [,4]
#> residual    1.0000000 4.0000000 16.000000 64.0000000
#> pseudohuber 0.4230609 0.7671455  0.828685  0.8330405
#> gaussian    0.6944444 2.7777778 11.111111 44.4444444

# The expected information is a quadrature here, and the family says so.
distributions7:::expected_hessian_exact(d)
#> [1] FALSE

# The quantile inverts the distribution function, and the generator
# inverts it at uniform variates, so a draw costs a root-find over a
# quadrature. That is what makes fitting this family dear.
q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
rbind(quantile = q, back = distrib_cdf(d, q, th))
#>                [,1] [,2]     [,3]
#> quantile -0.8235709  0.4 1.623571
#> back      0.2500000  0.5 0.750000
```
