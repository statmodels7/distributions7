# Expected Hessian of a Truncated Distribution

\\\mathbb{E}\[d\_{ij}\ell_T\] = -\mathrm{Cov}\_T(s_i, s_j)\\, the second
Bartlett identity under the truncated law.

## Usage

``` r
trunc_expected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  approx = c("bartlett", "integrate", "mc", "opg"),
  nsim = 10000,
  ...
)
```

## Arguments

- distrib:

  A truncated distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

## Value

A named list, one component per Hessian entry.

## Details

This still needs one quadrature per component even when the parent has
exact cdf derivatives: those give \\d^B Z\\ but cannot separate
\\\mathbb{E}\_T\[\ell^{(ij)}\]\\ from
\\\mathbb{E}\_T\[\ell^{(i)}\ell^{(j)}\]\\, and the covariance needs
both.
