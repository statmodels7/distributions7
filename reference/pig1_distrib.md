# Poisson-Inverse Gaussian Distribution Object

Builds a Poisson-inverse Gaussian distribution in its mean-dispersion
parametrization: mean \\\mu\\ and variance \\\mu + \sigma\mu^2\\, which
is gamlss's `PIG`. The family is the Poisson mixed over an inverse
Gaussian rate, an overdispersed count model with a **heavier upper tail
than the negative binomial at the same variance**.

## Usage

``` r
pig1_distrib(link_mu = log_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A `linkfunctions7` link object for the mean \\\mu\\, which must be
  strictly positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_sigma:

  A link object for the dispersion \\\sigma\\, also strictly positive.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class
[Pig1Distrib](https://statmodels7.github.io/distributions7/reference/Pig1Distrib.md),
inheriting from `discrete_distrib`. Its `params` are `c("mu", "sigma")`,
its `bounds` `c(0, Inf)`, and its `link_params` the two links given
here.

## The mass function, and how it is computed

The mass carries the modified Bessel function \\K\_{y-1/2}\\, which at
half-integer order is a finite sum. What the kernels evaluate is that
sum, on the log scale, after the prefactors have canceled; no Bessel
routine is called. With \\c = 1 + 2\sigma\mu\\ and \\\alpha =
\sqrt{c}/\sigma\\, \$\$\ell(y) = y\log\mu - \tfrac{y}{2}\log c +
\tfrac{1}{\sigma} + \psi(\alpha) - \log y!,\qquad \psi(\alpha) =
-\alpha + \log S_y(\alpha),\$\$ \$\$S_y(\alpha) = \sum\_{k=0}^{y-1}
\dfrac{\Gamma(y+k)}{\Gamma(k+1)\Gamma(y-k)}\\(2\alpha)^{-k}, \qquad S_0
= 1.\$\$ \\S_y\\ sums \\y\\ positive terms on the log scale, so nothing
cancels. The derivatives of \\\psi\\ in \\\alpha\\ are the weighted
rising-factorial moments of \\k\\ under those terms; everything else is
elementary. Every partial to fourth order is written out by hand, so no
order is differenced.

## One kernel per generic

There is one compiled kernel per quantity – `pig1_pdf_cpp`,
`pig1_gradient_cpp`, `pig1_hessian_cpp`, `pig1_deriv3_cpp`,
`pig1_deriv4_cpp` and their `pig2` twins – so asking for the mass does
not evaluate four orders of derivatives. Each recycles a scalar
parameter and tests the support itself, as every other compiled family's
kernel does: a \\y\\ that is negative, fractional or not finite reads
`-Inf` in the mass and `NaN` in a derivative.

The saving is smaller than the count of quantities suggests, because the
cost is \\S_y\\, which is \\O(y)\\ per observation and which every order
pays in full: measured over \\2\times10^4\\ observations the five
kernels run at 21.0, 21.3, 22.0, 22.7 and 23.3 milliseconds.

A second pair, `pig1_hd_jet_cpp` and `pig2_hd_jet_cpp`, carries a
bivariate jet truncated at total order four through the same expression.
It shares no algebra with the explicit route, so the tests compare the
two with no tolerance to hide behind, and it is not used in production;
`pig1_hd_cpp` and `pig2_hd_cpp` return all fifteen columns at once and
exist to be that comparison's other side.

## Notation

\\\mu\\ is the mean, \\\sigma\\ the dispersion, \\\alpha\\ the Bessel
argument, \\K\_\nu\\ the modified Bessel function of the second kind,
and \\\ell\\ the log-mass of one observation.

The **expected information** has no closed form and goes through the
summation strategies of
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).The
tail

Measured at \\\mu = 3\\, \\\sigma = 0.8\\, against a negative binomial
matched on the variance (\\\theta = 1/\sigma\\, both at 10.2): the mass
at \\y = 20\\ is \\6.8\times10^{-4}\\ against \\4.8\times10^{-4}\\, at
40 it is \\5.8\times10^{-6}\\ against \\5.4\times10^{-7}\\, and at 60
\\7.2\times10^{-8}\\ against \\5.6\times10^{-10}\\. The gap widens with
the count, and that widening is what a heavier tail means here.Which
parametrization to use

In this one \\\mu\\ and \\\sigma\\ are **not** orthogonal: measured at
the same setting, the mixed entry of the expected information summed
over the support is 7.39.
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
replaces \\\sigma\\ by the Bessel argument \\\alpha\\ and makes that
entry zero, at the cost of a second parameter with no moment reading.
Use this one to model the dispersion directly and that one to fit the
two parameters independently.Parameter domains

- \\\mu \in (0, \infty)\\

- \\\sigma \in (0, \infty)\\

## The distribution

\$\$P(Y=y) =
\sqrt{\frac{2\alpha}{\pi}}\\\frac{\mu^{y}e^{1/\sigma}}{(\alpha\sigma)^{y}\\y!}\\K\_{y-1/2}(\alpha),
\qquad \alpha = \sqrt{\frac{1}{\sigma^{2}} + \frac{2\mu}{\sigma}}\$\$ on
\\y \in \\0, 1, \dots\\\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \mu +
\sigma\mu^{2}\$\$

## References

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics* 54(3),
507–554.

Dean, C., Lawless, J. F., and Willmot, G. E. (1989). A mixed
Poisson-inverse-Gaussian regression model. *Canadian Journal of
Statistics* 17(2), 171–181.

## See also

[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
for the orthogonal parametrization,
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the other overdispersed count family,
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
for the limit at \\\sigma \to 0\\, the compiled kernels for the kernel,
and
[Pig1Distrib](https://statmodels7.github.io/distributions7/reference/Pig1Distrib.md)
for the class.

## Examples

``` r
d <- pig1_distrib()
th <- list(mu = 3, sigma = 0.8)

distrib_pdf(d, 0:5, th)
#> [1] 0.17197629 0.21422781 0.17775288 0.12895666 0.08968701 0.06196187
c(mean = mean(d, th), variance = variance(d, th),
  skewness = skewness(d, th))
#>      mean  variance  skewness 
#>  3.000000 10.200000  2.346499 

# Heavier in the tail than a negative binomial of the same variance.
nb <- negbin2_distrib()
nbth <- list(mu = 3, theta = 1 / 0.8)
c(variance_pig = variance(d, th), variance_negbin = variance(nb, nbth))
#>    variance_pig variance_negbin 
#>            10.2            10.2 
rbind(pig = distrib_pdf(d, c(20, 40, 60), th),
      negbin = distrib_pdf(nb, c(20, 40, 60), nbth))
#>                [,1]         [,2]         [,3]
#> pig    0.0006790258 5.753016e-06 7.233866e-08
#> negbin 0.0004803803 5.368144e-07 5.596723e-10

# The two parameters are not orthogonal, which pig2_distrib() repairs.
c(pig1 = sum(distrib_expected_hessian(d, 0:200, th,
                                      approx = "bartlett")$mu_sigma),
  pig2 = sum(distrib_expected_hessian(pig2_distrib(), 0:200,
                                      list(mu = 3, alpha = 3.010399),
                                      approx = "bartlett")$mu_alpha))
#>          pig1          pig2 
#>  7.392208e+00 -5.459445e-15 

# A fit recovers both parameters.
set.seed(63)
x <- distrib_rng(d, 4000, th)
coef(fit_distrib(d, x))
#>        mu     sigma 
#> 3.0102487 0.8074771 
```
