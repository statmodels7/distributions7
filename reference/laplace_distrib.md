# Laplace Distribution Object

Creates a distribution object for the Laplace (double-exponential)
distribution, parameterized by location (\\\mu\\) and scale
(\\\sigma\\). This is the reference example of a distribution whose
log-likelihood is not differentiable in a parameter (\\\mu\\); see
**Details** for how the package handles this.

## Usage

``` r
laplace_distrib(link_mu = identity_link(), link_sigma = log_link())
```

## Arguments

- link_mu:

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link function object for the scale parameter \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `LaplaceDistrib` (inheriting from
`continuous_distrib`).

## Details

The Laplace (double-exponential) distribution has location \\\mu\\ and
scale \\\sigma\\.

**Probability density function:** \$\$f(y; \mu, \sigma) =
\dfrac{1}{2\sigma} \exp\left(-\dfrac{\|y-\mu\|}{\sigma}\right)\$\$

**Cumulative distribution function:** \$\$F(q; \mu, \sigma) =
\begin{cases} \tfrac{1}{2}\exp\left(\tfrac{q-\mu}{\sigma}\right) & q \<
\mu \\ 1 - \tfrac{1}{2}\exp\left(-\tfrac{q-\mu}{\sigma}\right) & q \ge
\mu \end{cases}\$\$

**Quantile function:** \$\$Q(p; \mu, \sigma) = \mu -
\sigma\\\mathrm{sign}(p - \tfrac{1}{2})\\\log\left(1 - 2\left\|p -
\tfrac{1}{2}\right\|\right)\$\$

**Score** (defined almost everywhere): \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{\mathrm{sign}(y-\mu)}{\sigma}, \qquad
\dfrac{\partial \ell}{\partial \sigma} =
\dfrac{1}{\sigma}\left(\dfrac{\|y-\mu\|}{\sigma} - 1\right)\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = 0,
\quad \dfrac{\partial^2 \ell}{\partial \mu\\\partial \sigma} =
-\dfrac{\mathrm{sign}(y-\mu)}{\sigma^2}, \quad \dfrac{\partial^2
\ell}{\partial \sigma^2} = \dfrac{\sigma - 2\|y-\mu\|}{\sigma^3}\$\$

**Expected Hessian** (Fisher information from the score variance; see
below): \$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial
\mu^2}\right\] = -\dfrac{1}{\sigma^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{1}{\sigma^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu\\\partial \sigma}\right\] = 0\$\$

**Moments:** mean \\\mu\\, variance \\2\sigma^2\\, skewness 0, excess
kurtosis 3.

**Non-differentiability.** The density has a kink at \\y = \mu\\, so the
log-likelihood is not differentiable in \\\mu\\. The package marks this
via `params_smooth = c(mu = FALSE, sigma = TRUE)` and handles it as
follows:

- the **score**
  ([`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md))
  exists almost everywhere, equal to \\\mathrm{sign}(y-\mu)/\sigma\\;

- the **observed Hessian**
  ([`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md))
  has \\\partial^2 \ell / \partial \mu^2 = 0\\, so Newton-Raphson cannot
  update \\\mu\\;

- the **expected Hessian**
  ([`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md))
  is implemented in closed form from the variance of the score, giving
  the correct Fisher information \\1/\sigma^2\\ for \\\mu\\ and making
  Fisher scoring the appropriate estimation method. Because the closed
  form exists, the `approx` argument is ignored for this distribution.

**Parameter Domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

## Examples

``` r
d <- laplace_distrib()
d@params
#> [1] "mu"    "sigma"

theta <- list(mu = 0, sigma = 1)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.1839397 0.5000000 0.1839397
distrib_gradient(d, c(-1, 0, 1), theta)
#> $mu
#> [1] -1  0  1
#> 
#> $sigma
#> [1]  0 -1  0
#> 
```
