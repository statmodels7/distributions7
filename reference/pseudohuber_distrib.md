# Pseudo-Huber Distribution Object (Location-Scale Parameterization)

Creates a distribution object for the Pseudo-Huber distribution, whose
density corresponds to the exponential of the negative Pseudo-Huber
loss. It is a special case of the Generalized Hyperbolic distribution.

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

  A link function object for the location parameter \\\mu\\. Defaults to
  [`identity_link`](https://rdrr.io/pkg/linkfunctions7/man/identity_link.html).

- link_sigma:

  A link function object for the scale parameter \\\sigma\\. Defaults to
  [`log_link`](https://rdrr.io/pkg/linkfunctions7/man/log_link.html) to
  ensure positivity.

- link_nu:

  A link function object for the shape parameter \\\nu\\. Defaults to
  [`log_link`](https://rdrr.io/pkg/linkfunctions7/man/log_link.html) to
  ensure positivity.

## Value

An S7 object of class `PseudoHuberDistrib` (inheriting from
`continuous_distrib`) representing the Pseudo-Huber distribution.

## Details

The probability density function is: \$\$f(y; \mu, \sigma, \nu) =
\dfrac{1}{2 \sigma \sqrt{\nu} K_1(\sqrt{\nu})} \exp\left( - \sqrt{\nu +
\left(\dfrac{y-\mu}{\sigma}\right)^2} \right)\$\$ where \\K_1\\ is the
modified Bessel function of the second kind.

**Moments:**

- Expected value: \\\mathbb{E}(Y) = \mu\\

- Variance: \\\mathbb{V}(Y) = \sigma^2 \sqrt{\nu}\\
  \dfrac{K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})}\\

- Skewness: 0 (symmetric)

- Excess kurtosis: \\3 \dfrac{K_3(\sqrt{\nu})
  K_1(\sqrt{\nu})}{K_2(\sqrt{\nu})^2} - 3\\

**Score** (with \\r = y-\mu\\ and \\D = \sqrt{\nu + (r/\sigma)^2}\\;
\\K_1'\\ the derivative of \\K_1\\): \$\$\dfrac{\partial \ell}{\partial
\mu} = \dfrac{r}{\sigma^2 D}, \qquad \dfrac{\partial \ell}{\partial
\sigma} = \dfrac{1}{\sigma}\left(\dfrac{r^2}{\sigma^2 D} - 1\right)\$\$
\$\$\dfrac{\partial \ell}{\partial \nu} =
-\dfrac{1}{2}\left\[\dfrac{1}{\nu} + \dfrac{1}{D} +
\dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu}\\K_1(\sqrt{\nu})}\right\]\$\$ The
observed Hessian is available in closed form via
[`distrib_hessian.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md);
the expected Hessian has no closed form and is obtained by numerical
integration.

**Parameter Domains:**

- \\\mu \in (-\infty, +\infty)\\

- \\\sigma \in (0, +\infty)\\

- \\\nu \in (0, +\infty)\\

**Note:** The CDF, quantile function and RNG have no closed form and
rely on numerical integration / root-finding, so they are slower than
for the other distributions in the package. Response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are available in closed form; third- and fourth-order parameter
derivatives use the numerical fallback.

## See also

- [`distrib_pdf.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.PseudoHuberDistrib.md)
  for the density function.

- [`distrib_cdf.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuberDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.PseudoHuberDistrib.md)
  for the quantile function.

- [`distrib_rng.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.PseudoHuberDistrib.md)
  for random number generation.

- [`distrib_gradient.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuberDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.PseudoHuberDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PseudoHuberDistrib.md)
  for the expected Hessian.
