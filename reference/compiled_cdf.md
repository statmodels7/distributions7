# Distribution Function and Its Derivatives from the Compiled Kernels

Method bodies for
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md),
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
that read the compiled distribution functions of the package's scalar
registry, for the families whose distribution function has no closed
derivative in a shape parameter: the gamma (both parametrizations), the
chi-squared, the two generalized gammas, the two betas, the two von
Mises, the two Student t, the two pseudo-Huber and the skew t.

## Usage

``` r
compiled_cdf(distrib, q, theta, lower.tail = TRUE, log.p = FALSE, ...)

compiled_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)

compiled_hess_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)

compiled_cdf_args(distrib, q, theta)
```

## Arguments

- distrib:

  A distribution object of one of the families above.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameter values.

- lower.tail:

  Is the lower tail wanted? A single logical.

- ...:

  Unused.

- log, log.p:

  Are derivatives of (or values of) the log probability wanted? A single
  logical.

## Value

For `compiled_cdf()`, a numeric vector; for the two derivative bodies, a
named list of numeric vectors as
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
return.

`compiled_cdf_args()` returns the class name, the recycled quantiles and
the parameter matrix the kernels read.

## Details

A derivative in a shape parameter is the integral of the density's own
derivative, \$\$\partial_k F(q) = \int\_{-\infty}^{q} f\\ s_k \\ dy,
\qquad \partial\_{kl} F(q) = \int\_{-\infty}^{q} f\\(\ell\_{kl} + s_k
s_l)\\ dy,\$\$ with \\s\\ and \\\ell\\ the score and the Hessian of the
log-density. The integral is taken by a fixed rule (tanh-sinh panels
doubling away from the family's center, an exp-sinh tail beyond), in
\\\log y\\ on \\(0, \infty)\\ and in \\\mathrm{logit}(y)\\ on \\(0,
1)\\, so that a density that behaves like \\y^{k-1}\\ at an end of its
support is integrated as an exponential tail. Above the family's center
the integral runs over \\\[q, \infty)\\ with the opposite sign. For the
families that are location-scale in \\(\mu, \sigma)\\ only the shape
parameters are integrated: with \\z = (q - \mu)/\sigma\\ the location
and scale components are \\-f\\, \\-z f\\ and, at second order, \\-f
s\_\mu\\, \\f(-z^2 s\_\mu + 2z/\sigma)\\ and \\f(-z s\_\mu +
1/\sigma)\\, and the mixed components with a shape \\k\\ are \\-f s_k\\
and \\-z f s_k\\, all at \\q\\. The distribution function itself is
[`pt()`](https://rdrr.io/r/stats/TDist.html),
[`pgamma()`](https://rdrr.io/r/stats/GammaDist.html),
[`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) or
[`pbeta()`](https://rdrr.io/r/stats/Beta.html) where the family has one,
and the same integral of the density otherwise. The tail and the log
scale are applied as for every family, by the conversion of the
lower-tail derivatives.
