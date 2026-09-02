# Truncated Distribution Object

Restricts an existing distribution to \\\[L, U\]\\ and renormalizes it,
so that the probability the parent placed outside the interval is
redistributed inside it. Either endpoint may be omitted, giving
one-sided truncation, and at least one must be given. Discrete and
continuous parents are both accepted.

BOTH ENDPOINTS ARE INCLUDED, which for a discrete parent is the
difference between `truncated(poisson_distrib(), lower = 1)`, the
zero-truncated Poisson on \\\\1, 2, \dots\\\\, and truncating above one.

## Usage

``` r
truncated(distrib, lower = NULL, upper = NULL)
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` or `continuous_distrib`,
  including one that is already truncated.

- lower, upper:

  The truncation points \\L\\ and \\U\\, each a single number or `NULL`
  for no truncation on that side. At least one must be supplied, and for
  a discrete parent a finite one must be a whole number. Both are
  included in the resulting support.

## Value

An S7 object of class
[TruncatedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedDiscreteDistrib.md)
or
[TruncatedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedContinuousDistrib.md),
matching the parent's branch. Its `params`, `params_bounds` and
`link_params` are the parent's unchanged; its `bounds` are \\\[L, U\]\\,
and `distrib_name` is `"truncated "` followed by the parent's name and
the interval.

## The three functions

Write \\Z(\theta) = P(L \le Y \le U)\\ for the retained mass. Then
\$\$f_T(y;\theta) = \frac{f(y;\theta)}{Z(\theta)} \quad (L \le y \le U),
\qquad F_T(q) = \frac{F(q;\theta) - F(L^-;\theta)}{Z(\theta)}, \qquad
Q_T(p) = Q\\F(L^-;\theta) + pZ(\theta)\\,\$\$ with \\F(L^-) = F(L)\\ for
a continuous parent and \\F(L) - f(L)\\ for a discrete one, the lower
endpoint being kept. The quantile function inverts \\F_T\\ through the
parent's own, so drawing from the truncated law is an exact one-pass
inverse transform however small \\Z\\ is.

## Truncation adds no parameter

The endpoints are known constants, like a binomial's `size`, so the
result carries exactly the parent's parameters, domains and links. What
it adds is the \\\theta\\-dependent normalizing constant, which
contributes to every derivative of \\\ell_T = \ell - \log Z\\. Writing
\$\$m_i = \mathbb{E}\_T\[s_i\], \qquad M\_{ij} =
\mathbb{E}\_T\[H\_{ij} + s_i s_j\],\$\$ with both expectations taken
under the TRUNCATED law, \$\$\frac{\partial \ell_T}{\partial\theta_i} =
s_i(y) - m_i, \qquad
\frac{\partial^{2}\ell_T}{\partial\theta_i\partial\theta_j} =
H\_{ij}(y) - M\_{ij} + m_i m_j,\$\$
\$\$\mathbb{E}\left\[\frac{\partial^{2}\ell_T}
{\partial\theta_i\partial\theta_j}\right\] = -\mathrm{Cov}\_T(s_i,
s_j).\$\$ The score is the parent's recentered at its truncated mean,
and the information is the covariance of that score, which is the second
Bartlett identity for the truncated model and is used here as a
consistency check on the two expressions above.

## What this costs

\\m_i\\ and \\M\_{ij}\\ have no closed form for a general parent, and
are obtained by quadrature or summation through
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).
Where the parent has exact cdf derivatives they come instead from two
calls on it, measured at 1.4 ms against 4.9 ms for a truncated
gaussian's Hessian; see
[`has_exact_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/has_exact_cdf_deriv.md)
for why that route is gated on accuracy. Third and fourth derivatives
are closed form as well, assembled by the partition sums of
`wrapper_derivatives.R` from the same ratios taken at higher blocks, and
they take the same two routes: the parent's cdf derivative of that order
where the gate admits it, and one memoized
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
per distinct block otherwise.

## What the constructor rejects

- Both endpoints `NULL`. There is nothing to do, and returning the
  parent silently would hide the missing argument.

- `lower >= upper`.

- An endpoint removing no mass, such as
  `truncated(gamma2_distrib(), lower = -2)`: the gamma lives on
  \\(0,\infty)\\, so the result would be the gamma itself.

- A fractional endpoint on a discrete parent, which is ambiguous between
  the support points either side of it.

- A discrete interval leaving too few support points to identify the
  parameters: \\k\\ points carry \\k-1\\ free probabilities, so
  `n_params + 1` points are required.

- A parent modeling a probability of zero, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
  or
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
  when the interval removes \\0\\. Truncating zero away cancels that
  parameter out of the likelihood, leaving an identically zero score.
  Truncating elsewhere, as in
  `truncated(zero_adjusted(gamma2_distrib()), upper = 5)`, is legitimate
  and the point mass is carried through
  [`distrib_atoms.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.TruncatedContinuousDistrib.md).

Truncating an already truncated distribution is allowed and is COLLAPSED
into one object over the intersection of the two intervals. Nesting
would be correct and would pay the quadrature cost twice.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\,
\\F\\ and \\Q\\ are the parent's density, distribution and quantile
functions; \\s_i\\ and \\H\_{ij}\\ are the parent's score and observed
Hessian; \\\ell_T\\ is the truncated log-likelihood; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the wrappers that DO add a parameter,
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
for a change of variable,
[TruncatedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedContinuousDistrib.md)
and
[TruncatedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/TruncatedDiscreteDistrib.md)
for the two classes, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate the result.

## Examples

``` r
# The zero-truncated Poisson: the lower endpoint is kept.
ztp <- truncated(poisson_distrib(), lower = 1)
distrib_pdf(ztp, 0:4, list(mu = 2))
#> [1] 0.0000000 0.3130353 0.3130353 0.2086902 0.1043451
dpois(1:4, 2) / (1 - dpois(0, 2))
#> [1] 0.3130353 0.3130353 0.2086902 0.1043451
c(fitted_mean = mean(ztp, list(mu = 2)), theory = 2 / (1 - exp(-2)))
#> fitted_mean      theory 
#>    2.313035    2.313035 

# A gaussian restricted to an interval. Truncation adds no parameter.
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)
identical(tn@params, gaussian1_distrib()@params)
#> [1] TRUE
c(mean = mean(tn, theta), variance = variance(tn, theta))
#>      mean  variance 
#> 0.4159512 0.6028535 

# The score is the parent's, recentered so that it has mean zero.
set.seed(1)
y <- distrib_rng(tn, 20000, theta)
round(vapply(distrib_gradient(tn, y, theta), mean, numeric(1)), 3)
#>    mu sigma 
#> 0.001 0.004 

# Nesting collapses to the intersection rather than paying twice.
t2 <- truncated(truncated(gaussian1_distrib(), lower = -1), upper = 2)
c(lower = t2@lower, upper = t2@upper)
#> lower upper 
#>    -1     2 

# A zero wrapper may be truncated elsewhere, and its atom is carried through.
tz <- truncated(zero_adjusted(gamma2_distrib()), upper = 5)
distrib_atoms(tz, list(mu = 2, sigma2 = 1, za = 0.3))
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3021864
#> 

# Three refusals, each naming the condition that failed.
try(truncated(gamma2_distrib(), lower = -2))
#> Error : truncated() was given lower = -2, which is at or below the lower bound of the
#>   support of 'gamma2' (0). Truncating there removes no probability mass and the
#>   result would be the parent distribution. Choose a point strictly inside the
#>   support, or omit 'lower'.
try(truncated(poisson_distrib(), lower = 1.5))
#> Error : 'lower' = 1.5 is not a point of the support. A discrete distribution is supported on
#>   the integers, so a non-integer truncation point is ambiguous: use 1 or 2.
try(truncated(zero_inflated(poisson_distrib()), lower = 1))
#> Error : truncated() cannot remove 0 from the support of 'zero-inflated poisson', which models the
#>   probability of a zero. The truncation constant cancels that parameter out of
#>   the likelihood, leaving an identically zero score. Truncate elsewhere, or
#>   truncate the underlying distribution before wrapping it.
```
