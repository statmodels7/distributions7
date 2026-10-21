# Distribution Function and Its Derivatives from the Compiled Kernels

Provides the method bodies for
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

compiled_log_cdf(distrib, a, out, lower.tail)

log_tail_integral(distrib, q, th, lower.tail)

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

  A single logical; `TRUE` selects the lower tail and `FALSE` the upper
  tail.

- ...:

  Unused.

- a:

  The value of `compiled_cdf_args()`.

- out:

  The probabilities of the requested tail, from the compiled kernel.

- th:

  A named list of scalar parameter values.

- log, log.p:

  A single logical; if `TRUE`, derivatives of (or values of) the log
  probability are returned.

## Value

For `compiled_cdf()`, a numeric vector; for the two derivative bodies, a
named list of numeric vectors as
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
return.

`compiled_far_tail()` returns the logical vector of the quantiles whose
requested tail is below \\e^{-5}\\, inside the support.

`compiled_log_cdf()` returns the logarithm of `out`, formed as
\\\log(1 - P\_{\mathrm{other}})\\ through
[`log1p()`](https://rdrr.io/r/base/Log.html) where `out` is above one
half, and as \\\ell(q) + \log\int\_{\mathrm{tail}} e^{\ell(y) -
\ell(q)}\\dy\\ where `out` is below \\e^{-5}\\; there
[`quad_log_tail_prob()`](https://statmodels7.github.io/distributions7/reference/quad_log_tail_derivs.md)
computes the integral on the scale \\1/\|\partial_y \ell(q)\|\\.

`log_tail_integral()` returns the logarithm of the probability of the
tail beyond `q`, \\\ell(q) + \log\int_0^\infty e^{\ell(q \mp u) -
\ell(q)}\\du\\, by
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) on a
scale adjusted until one step changes \\\ell\\ by between 1/2 and 50,
and \\\ell(q) - \log\|\ell'(q)\|\\ where such a step cannot be resolved.

`compiled_cdf_args()` returns the class name, the recycled quantiles and
the parameter matrix that the kernels read.

## Details

A derivative in a shape parameter is the integral of the density's own
derivative, \$\$\partial_k F(q) = \int\_{-\infty}^{q} f\\ s_k \\ dy,
\qquad \partial\_{kl} F(q) = \int\_{-\infty}^{q} f\\(\ell\_{kl} + s_k
s_l)\\ dy,\$\$ with \\s\\ and \\\ell\\ the score and the Hessian of the
log-density. The integral is taken by a fixed rule (tanh-sinh panels
doubling away from the family's center, an exp-sinh tail beyond), in
\\\log y\\ on \\(0, \infty)\\ and in \\\mathrm{logit}(y)\\ on \\(0,
1)\\, so that a density that behaves like \\y^{k-1}\\ at an end of its
support is integrated as an exponential tail; in \\\log y\\ the panels
are at most one unit wide, and they restart at the scale \\\theta\\
where the density's fast decay begins when that point lies more than one
unit above the center. For a small \\k\\ part of the mass lies below the
smallest representable double, so on \\(0, \infty)\\ the part below
\\y_0 = \min(q, 10^{-290})\\ is integrated in closed form wherever it
holds at least \\10^{-17}\\ of \\F(q)\\: there the log-density and its
parameter derivatives are linear in \\L = \log y\\, so they are read at
\\10^{-290}\\ and 15 units of \\L\\ below it and extended to \\L_0 =
\log y_0\\, and \\\int\_{-\infty}^{L_0} e^{c + k u} g(u)\\ du = f(y_0)\\
y_0 \sum_j (-1)^j g^{(j)}(L_0)/k^{j+1}\\
([`compiled_cdf_head()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf_head.md)
at orders three and four). When the part of the support above \\q\\
holds less mass than the part below, the integral runs over the upper
part with the opposite sign, so that a derivative that is small beside
one is not a sum of terms of order one; the mass is the family's closed
distribution function where it has one, and the integral of the density
otherwise. The two beta families always integrate the lower part when
their second shape is below one. For the families that are
location-scale in \\(\mu, \sigma)\\ only the shape parameters are
integrated: with \\z = (q - \mu)/\sigma\\ the location and scale
components are \\-f\\, \\-z f\\ and, at second order, \\-f s\_\mu\\,
\\f(-z^2 s\_\mu + 2z/\sigma)\\ and \\f(-z s\_\mu + 1/\sigma)\\, and the
mixed components with a shape \\k\\ are \\-f s_k\\ and \\-z f s_k\\, all
at \\q\\. The distribution function itself is
[`pt()`](https://rdrr.io/r/stats/TDist.html),
[`pgamma()`](https://rdrr.io/r/stats/GammaDist.html),
[`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) or
[`pbeta()`](https://rdrr.io/r/stats/Beta.html) where the family has one,
and the same integral of the density otherwise; with `log.p = TRUE` its
logarithm is formed by `compiled_log_cdf()`, so that it stays finite in
either tail. The lower-tail derivatives are converted to the requested
tail and log scale, as for every family.
