# Beta Distribution Object (Mean-Precision Parameterization)

Creates a distribution object for the Beta distribution parameterized by
mean (\\\mu\\) and precision (\\\phi\\).

## Usage

``` r
beta_distrib(link_mu = logit_link(), link_phi = log_link())
```

## Arguments

- link_mu:

  A link function object for the mean parameter \\\mu\\. Defaults to
  [`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html)
  to ensure the parameter stays within (0, 1).

- link_phi:

  A link function object for the precision parameter \\\phi\\. Defaults
  to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `BetaDistrib` (inheriting from
`continuous_distrib`) representing the Beta distribution.

## Details

The Beta distribution is given a mean/precision parameterization:
\\\mu\\ is the mean and \\\phi\\ a precision parameter. The standard
shapes are \$\$\alpha = \mu\phi, \qquad \beta = (1-\mu)\phi\$\$

**Probability density function:** \$\$f(y; \mu, \phi) =
\dfrac{\Gamma(\phi)}{\Gamma(\alpha)\Gamma(\beta)}\\ y^{\alpha-1}
(1-y)^{\beta-1}, \quad 0 \< y \< 1\$\$

**Cumulative distribution function** (\\I_q\\ the regularized incomplete
beta function): \$\$F(q; \mu, \phi) = I_q(\alpha, \beta)\$\$

**Quantile function:** no closed form; the numerical inverse of the CDF.

**Score** (\\\psi\\ the digamma function): \$\$\dfrac{\partial
\ell}{\partial \mu} = \phi\left\[\log\left(\dfrac{y}{1-y}\right) -
\psi(\alpha) + \psi(\beta)\right\]\$\$ \$\$\dfrac{\partial
\ell}{\partial \phi} = \psi(\phi) - \mu\psi(\alpha) -
(1-\mu)\psi(\beta) + \mu\log y + (1-\mu)\log(1-y)\$\$

**Expected Hessian** (\\\psi_1\\ the trigamma function; the observed and
expected Hessians coincide, as they do not depend on \\y\\):
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
-\phi^2\left\[\psi_1(\alpha) + \psi_1(\beta)\right\], \qquad
\dfrac{\partial^2 \ell}{\partial \mu\\\partial \phi} =
-\phi\left\[\mu\psi_1(\alpha) - (1-\mu)\psi_1(\beta)\right\]\$\$
\$\$\dfrac{\partial^2 \ell}{\partial \phi^2} = \psi_1(\phi) -
\mu^2\psi_1(\alpha) - (1-\mu)^2\psi_1(\beta)\$\$

**Moments:** mean \\\mu\\, variance \\\mu(1-\mu)/(\phi+1)\\, skewness
\\\dfrac{2(1-2\mu)\sqrt{\phi+1}}{(\phi+2)\sqrt{\mu(1-\mu)}}\\.

**Parameter domains:**

- \\\mu \in (0, 1)\\

- \\\phi \in (0, +\infty)\\

Analytical third- and fourth-order derivatives
([`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
and response derivatives
([`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md))
are also available.

## See also

- [`distrib_pdf.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaDistrib.md)
  for the probability density function.

- [`distrib_cdf.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BetaDistrib.md)
  for the cumulative distribution function.

- [`distrib_quantile.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BetaDistrib.md)
  for the quantile function.

- [`distrib_rng.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BetaDistrib.md)
  for random number generation.

- [`distrib_gradient.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaDistrib.md)
  for the analytical gradient.

- [`distrib_hessian.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaDistrib.md)
  for the analytical observed Hessian.

- [`distrib_expected_hessian.BetaDistrib`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaDistrib.md)
  for the analytical expected Hessian.

## Examples

``` r
d <- beta_distrib()
d@params
#> [1] "mu"  "phi"

theta <- list(mu = 0.4, phi = 5)
distrib_pdf(d, c(0.2, 0.5, 0.8), theta)
#> [1] 1.536 1.500 0.384
distrib_gradient(d, c(0.2, 0.5, 0.8), theta)
#> $mu
#> [1] -4.431472  2.500000  9.431472
#> 
#> $phi
#> [1]  0.005672038  0.090186153 -0.271586835
#> 
```
