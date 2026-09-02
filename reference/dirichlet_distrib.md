# Dirichlet Distribution, Mean Vector and Concentration

Builds the distribution object for the Dirichlet family on the open
simplex of \\p\\ coordinates, parametrized by a mean vector \\\mu\\ and
a concentration \\\phi \> 0\\. The returned object carries closed-form
derivatives of the log-density to fourth order and a **closed-form
expected information**, which is unusual for a family with no location
and scale to separate.

The mean is carried on a `parameters7` simplex and flattened into scalar
parameters, so every generic of the package indexes it as it does any
other family.

## Usage

``` r
dirichlet_distrib(
  n_dim,
  mean = parameters7::simplex(n_dim),
  link_phi = log_link()
)
```

## Arguments

- n_dim:

  The number of coordinates \\p\\, a single integer of at least 2.
  Anything else signals an error naming the argument.

- mean:

  A `parameters7` parameter producing \\p\\ coordinates that sum to one,
  normally a
  [`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html).
  Defaults to `parameters7::simplex(n_dim)`. An object that is not a
  `parameters7` parameter, or that produces a different number of
  coordinates, signals an error naming both counts.

- link_phi:

  A `link` object from `linkfunctions7` for the concentration \\\phi\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive. The mean's free values are unconstrained already and
  carry the identity.

## Value

An S7 object of class `DirichletDistrib`, inheriting from
`multivariate_distrib`, with `param` the simplex given here,
`distrib_name` `"dirichlet [pd, mean=<chart>]"`, `dimension`
`"multivariate"`, `n_dim` \\p\\, `bounds` `c(0, 1)`, `params` the
simplex's free names prefixed by `mean_` followed by `"phi"`, `n_params`
\\p\\, and `link_params` the identity for each mean coordinate and
`link_phi` for the concentration.

## The parametrization

Writing \\\alpha = \phi\mu\\ for the shapes, the density on the simplex
\\\\y_j \> 0, \sum_j y_j = 1\\\\ is \$\$f(y) =
\frac{\Gamma(\phi)}{\prod\_{j=1}^{p}\Gamma(\alpha_j)} \prod\_{j=1}^{p}
y_j^{\alpha_j - 1},\$\$ with \$\$\mathbb{E}\[Y_j\] = \mu_j, \qquad
\operatorname{Cov}(Y_i, Y_j) = \frac{\delta\_{ij}\mu_i -
\mu_i\mu_j}{\phi + 1}.\$\$ So \\\phi\\ is a **precision**: the larger it
is, the tighter the draws sit about the mean. The covariance is
singular, the coordinates summing to one, and every off-diagonal entry
is negative.

The parametrization follows the design of the multivariate gaussian's.
The constrained object, here a point of the simplex, is carried by a
`parameters7` parameter and **flattened into scalars** with identity
links, so
[`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md),
the derivative names, the link scale and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
need no special case. The parameter names are the simplex's own free
names prefixed by `mean_`, followed by `phi`.

## A multivariate family that is not elliptical

This is the first family of the package with no location and scale to
separate: the support is a set of dimension \\p-1\\, not a Euclidean
space, and the covariance is singular by construction. It is therefore
the real test of the multivariate layer, and two of its methods exist to
answer that test.
[`mv_marginal.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.DirichletDistrib.md)
returns a beta object rather than an error, and
[`mv_reference_draw.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.DirichletDistrib.md)
replaces the base class's gaussian proposal with the uniform on the
simplex, without which
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
would estimate an integral that is 1 at about `2e-08`.

## Why the expected information is closed form

Two identities do the work. Differentiating \\\sum_j \mu_j = 1\\ once
shows that the columns of \\A = \partial\mu/\partial\eta\\ sum to zero,
and twice that every second-derivative vector of the simplex does; and
\\\mathbb{E}\[\log y_j\] = \psi(\alpha_j) - \psi(\phi)\\ makes
\\\mathbb{E}\[g_j\] = -\psi(\phi)\\ the same constant for every \\j\\.
Every term of the observed Hessian that carries the data is that
constant times one of the two zero sums, and vanishes.

## What is refused

The distribution function and the quantile are refused by
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md),
as for every family of that class: a distribution function on
\\\mathbb{R}^p\\ is an orthant probability and a quantile needs an
ordering. The response derivatives are refused too. A marginal over
several coordinates at once is refused with the reason.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale, over the simplex's free
values and \\\log\phi\\. No estimate is closed form.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean vector,
\\\phi \> 0\\ the concentration, \\\alpha = \phi\mu\\ the shapes, \\p\\
the number of coordinates, \\\eta\\ the free vector of the simplex chart
and \\\psi\\ the digamma function.

## References

Kotz, S., Balakrishnan, N. and Johnson, N. L. (2000). *Continuous
Multivariate Distributions*, Volume 1, 2nd edition, Chapter 49. Wiley,
New York.

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall, London.

## See also

[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for a coordinate's marginal and the two-coordinate case seen on the
line;
[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
for the discrete family on the same simplex;
[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
for the elliptical multivariate family;
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
for the chart the mean rides;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[DirichletDistrib](https://statmodels7.github.io/distributions7/reference/DirichletDistrib.md)
for the class.

## Examples

``` r
d <- dirichlet_distrib(3)
d
#> Distribution: Dirichlet [3d, Mean=simplex]
#> Type:         Continuous, 3-dimensional
#> Dimensions:   multivariate
#> 
#> Parameters:
#>   mean_alr1 (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   mean_alr2 (mean)               | Link: identity   | Domain: (-Inf, Inf)
#>   phi       (concentration)      | Link: log        | Domain: (0, Inf)

# Two free mean values and a concentration, all on the parameter scale.
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
mv_location(d, th)
#> [1] 0.4260125 0.2583897 0.3155978
round(mv_sigma(d, th), 6)
#>           [,1]      [,2]      [,3]
#> [1,]  0.018810 -0.008467 -0.010342
#> [2,] -0.008467  0.014740 -0.006273
#> [3,] -0.010342 -0.006273  0.016615

# The covariance is singular: rank p - 1, with the ones vector in the null
# space.
c(rank = qr(mv_sigma(d, th))$rank, dim = 3)
#> rank  dim 
#>    2    3 

# The marginals are beta with the same concentration, so a panel of a pairs
# plot is a real distribution.
mv_marginal(d, th, which = 1)$theta
#> $mu
#> [1] 0.4260125
#> 
#> $phi
#> [1] 12
#> 

# The density and the sample agree on the mean and the coordinate variance.
set.seed(3)
Z <- distrib_rng(d, 3e5, th)
rbind(sample = c(mean(Z[, 1]), var(Z[, 1])),
      theoretical = c(mv_location(d, th)[1], mv_sigma(d, th)[1, 1]))
#>                  [,1]       [,2]
#> sample      0.4260478 0.01883675
#> theoretical 0.4260125 0.01880968

# Fitting recovers the parameters.
set.seed(9)
coef(fit_distrib(d, distrib_rng(d, 800, th)))
#>  mean_alr1  mean_alr2        phi 
#>  0.3188722 -0.2176956 12.4357208 
```
