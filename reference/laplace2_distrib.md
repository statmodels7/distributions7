# Laplace Distribution in Location and Rate

Creates a Laplace (double-exponential) distribution object parametrized
by location (\\\mu\\) and **rate** (\\\lambda\\).

## Usage

``` r
laplace2_distrib(link_mu = identity_link(), link_lambda = log_link())
```

## Arguments

- link_mu:

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_lambda:

  A link function object for the rate parameter \\\lambda\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `Laplace2Distrib` (inheriting from
`continuous_distrib`).

## Details

This is the same law as
[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
in different coordinates: \\\lambda\\ here is \\1/\sigma\\ there. The
rate form is the one a penalty consumes, because the negative
log-density at fixed \\\mu = 0\\ is \\\lambda \lvert y \rvert\\ up to a
constant: the lasso penalty is linear in \\\lambda\\, and its
derivatives in \\\lambda\\ beyond the first carry no data at all.

**Probability density function:** \$\$f(y; \mu, \lambda) =
\dfrac{\lambda}{2} \exp\left(-\lambda\|y-\mu\|\right)\$\$

**Score** (defined almost everywhere): \$\$\dfrac{\partial
\ell}{\partial \mu} = \lambda\\\mathrm{sign}(y-\mu), \qquad
\dfrac{\partial \ell}{\partial \lambda} = \dfrac{1}{\lambda} -
\|y-\mu\|\$\$

**Observed Hessian** (almost everywhere): \$\$\dfrac{\partial^2
\ell}{\partial \mu^2} = 0, \quad \dfrac{\partial^2 \ell}{\partial
\mu\\\partial \lambda} = \mathrm{sign}(y-\mu), \quad \dfrac{\partial^2
\ell}{\partial \lambda^2} = -\dfrac{1}{\lambda^2}\$\$

**Expected Hessian** (Fisher information from the score variance):
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\lambda^2, \quad \mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial
\lambda^2}\right\] = -\dfrac{1}{\lambda^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu\\\partial
\lambda}\right\] = 0\$\$

**Moments:** mean \\\mu\\, variance \\2/\lambda^2\\, skewness 0, excess
kurtosis 3.

The kink at \\y = \mu\\ and its consequences are those of
[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md):
`params_smooth = c(mu = FALSE, lambda = TRUE)`, the observed Hessian
cannot update \\\mu\\, and the expected Hessian comes from the variance
of the score.

**Parameter Domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\lambda \in (0, +\infty)\\

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)

## Examples

``` r
d <- laplace2_distrib()
theta <- list(mu = 0, lambda = 2)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.1353353 1.0000000 0.1353353

# the same law as laplace_distrib with sigma = 1/2
distrib_pdf(laplace_distrib(), 0.7, list(mu = 0, sigma = 0.5))
#> [1] 0.246597
```
