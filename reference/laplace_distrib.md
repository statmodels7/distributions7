# Laplace Distribution Object

Creates a distribution object for the Laplace (double-exponential)
distribution, parameterized by location (\\\mu\\) and scale (\\b\\).
This is the reference example of a distribution whose log-likelihood is
not differentiable in a parameter (\\\mu\\); see **Details** for how the
package handles this.

## Usage

``` r
laplace_distrib(link_mu = identity_link(), link_b = log_link())
```

## Arguments

- link_mu:

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_b:

  A link function object for the scale parameter \\b\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `LaplaceDistrib` (inheriting from
`continuous_distrib`).

## Details

The Laplace (double-exponential) distribution has location \\\mu\\ and
scale \\b\\.

**Probability density function:** \$\$f(y; \mu, b) = \dfrac{1}{2b}
\exp\left(-\dfrac{\|y-\mu\|}{b}\right)\$\$

**Cumulative distribution function:** \$\$F(q; \mu, b) = \begin{cases}
\tfrac{1}{2}\exp\left(\tfrac{q-\mu}{b}\right) & q \< \mu \\ 1 -
\tfrac{1}{2}\exp\left(-\tfrac{q-\mu}{b}\right) & q \ge \mu
\end{cases}\$\$

**Quantile function:** \$\$Q(p; \mu, b) = \mu - b\\\mathrm{sign}(p -
\tfrac{1}{2})\\\log\left(1 - 2\left\|p - \tfrac{1}{2}\right\|\right)\$\$

**Score** (defined almost everywhere): \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{\mathrm{sign}(y-\mu)}{b}, \qquad
\dfrac{\partial \ell}{\partial b} =
\dfrac{1}{b}\left(\dfrac{\|y-\mu\|}{b} - 1\right)\$\$

**Observed Hessian:** \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = 0,
\quad \dfrac{\partial^2 \ell}{\partial \mu\\\partial b} =
-\dfrac{\mathrm{sign}(y-\mu)}{b^2}, \quad \dfrac{\partial^2
\ell}{\partial b^2} = \dfrac{b - 2\|y-\mu\|}{b^3}\$\$

**Expected Hessian** (Fisher information from the score variance; see
below): \$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial
\mu^2}\right\] = -\dfrac{1}{b^2}, \quad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial b^2}\right\] =
-\dfrac{1}{b^2}, \quad \mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial
\mu\\\partial b}\right\] = 0\$\$

**Moments:** mean \\\mu\\, variance \\2b^2\\, skewness 0, excess
kurtosis 3.

**Non-differentiability.** The density has a kink at \\y = \mu\\, so the
log-likelihood is not differentiable in \\\mu\\. The package marks this
via `params_smooth = c(mu = FALSE, b = TRUE)` and handles it as follows:

- the **score**
  ([`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md))
  exists almost everywhere, equal to \\\mathrm{sign}(y-\mu)/b\\;

- the **observed Hessian**
  ([`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md))
  has \\\partial^2 \ell / \partial \mu^2 = 0\\, so Newton-Raphson cannot
  update \\\mu\\;

- the **expected Hessian**
  ([`distrib_expected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md))
  is defined from the score variance (first Bartlett identity), giving
  the correct Fisher information \\1/b^2\\ for \\\mu\\ and making Fisher
  scoring the appropriate estimation method.

**Parameter Domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\b \in (0, +\infty)\\
