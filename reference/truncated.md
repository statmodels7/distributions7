# Truncated Distribution Object

Restricts an existing distribution to \\\[\ell, u\]\\ and renormalises
it, so that all the probability mass the parent placed outside the
interval is redistributed inside it. Either endpoint may be omitted,
giving one-sided truncation; at least one must be given.

Works on discrete and continuous parents alike. **Both endpoints are
included**, which for a discrete parent is the difference between
`truncated(poisson_distrib(), lower = 1)` — the zero-truncated Poisson,
supported on \\\\1, 2, \dots\\\\ — and truncating above 1.

## Usage

``` r
truncated(distrib, lower = NULL, upper = NULL)
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` or `continuous_distrib`.

- lower, upper:

  The truncation points. Each may be `NULL` (no truncation on that
  side), and at least one must be supplied. For a discrete parent both
  must be whole numbers.

## Value

An S7 object of class `TruncatedDiscreteDistrib` or
`TruncatedContinuousDistrib`.

## Details

Write \\Z(\theta) = P(\ell \le Y \le u)\\ for the retained mass. Then
\$\$f_T(y;\theta) = \frac{f(y;\theta)}{Z(\theta)}\\ \\ (\ell \le y \le
u), \qquad F_T(q) = \frac{F(q;\theta) - F(\ell^-;\theta)}{Z(\theta)},
\qquad Q_T(p) = Q\\\left(F(\ell^-;\theta) + pZ(\theta)\right),\$\$ with
\\F(\ell^-) = F(\ell)\\ for a continuous parent and \\F(\ell) -
f(\ell)\\ for a discrete one, since the lower endpoint is kept.

**Truncation adds no parameter.** The endpoints are known constants,
like a binomial's `size`, so the truncated distribution has exactly the
parent's parameters, domains and links. What it does add is a
\\\theta\\-dependent normalising constant, and that constant contributes
to every derivative of the log-likelihood \\\ell_T = \ell - \log Z\\.
Writing \$\$m_i = \mathbb{E}\_T\[s_i\], \qquad M\_{ij} =
\mathbb{E}\_T\[H\_{ij} + s_is_j\],\$\$ where the expectations are taken
under the *truncated* distribution, \$\$\frac{\partial
\ell_T}{\partial\theta_i} = s_i(y) - m_i, \qquad
\frac{\partial^{2}\ell_T}{\partial\theta_i\partial\theta_j} =
H\_{ij}(y) - M\_{ij} + m_im_j,\$\$
\$\$\mathbb{E}\left\[\frac{\partial^{2}\ell_T}{\partial\theta_i\partial\theta_j}\right\]
= -\operatorname{Cov}\_T(s_i, s_j).\$\$ The score is simply the parent's
score *recentred* at its truncated mean, and the information is the
covariance of that score — which is the second Bartlett identity for the
truncated model, and is used as a consistency check rather than derived
separately.

**What this costs.** \\m_i\\ and \\M\_{ij}\\ have no closed form for a
general parent, and are obtained by quadrature (continuous) or summation
(discrete) through
[`expectation`](https://statmodels7.github.io/distributions7/reference/expectation.md).
Derivatives of a truncated distribution are therefore substantially more
expensive than the parent's, and third and fourth derivatives fall back
to finite differences of the analytical Hessian.

**What the constructor refuses.**

- Both endpoints `NULL`: nothing to do, and silently returning the
  parent would hide the mistake.

- `lower >= upper`.

- A truncation point that removes no mass, such as
  `truncated(gamma2_distrib(), lower = -2)`: the Gamma is supported on
  \\(0,\infty)\\, so the result would be the Gamma itself.

- A non-integer endpoint for a discrete parent, which is ambiguous.

- A discrete truncation leaving too few support points to identify the
  parameters: \\k\\ points carry \\k-1\\ free probabilities, so
  `n_params + 1` points are needed.

- A parent that models a probability of zero —
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
  or
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
  — when the truncation removes \\0\\ from the support. Truncating zero
  away cancels that parameter out of the likelihood entirely, leaving an
  identically zero score. Truncating elsewhere, as in
  `truncated(zero_adjusted(gamma2_distrib()), upper = 5)`, is fine and
  the point mass is carried through
  [`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md).

Truncating an already truncated distribution is allowed and is collapsed
into a single object over the intersection of the two intervals, rather
than nested.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md),
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)

## Examples

``` r
# The zero-truncated Poisson
ztp <- truncated(poisson_distrib(), lower = 1)
distrib_pdf(ztp, 0:4, list(mu = 2))
#> [1] 0.0000000 0.3130353 0.3130353 0.2086902 0.1043451
dpois(1:4, 2) / (1 - dpois(0, 2))
#> [1] 0.3130353 0.3130353 0.2086902 0.1043451

# A Gaussian restricted to an interval
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
mean(tn, list(mu = 0, sigma = 1))
#> [1] 0.2296372

# A truncation point that removes nothing is refused
try(truncated(gamma2_distrib(), lower = -2))
#> Error : truncated() was given lower = -2, which is at or below the lower bound of the
#>   support of 'gamma2' (0). Truncating there removes no probability mass and the
#>   result would be the parent distribution. Choose a point strictly inside the
#>   support, or omit 'lower'.
```
