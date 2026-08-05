# Dirichlet Distribution Object

Creates a distribution object for the Dirichlet distribution on the
simplex, parametrised by a mean vector and a concentration \\\phi\\.

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

  The number of coordinates \\p\\.

- mean:

  A parameters7
  [`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.html)
  of the same dimension, carrying the mean. Defaults to
  `parameters7::simplex(n_dim)`.

- link_phi:

  A link function object for \\\phi\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `DirichletDistrib`.

## Details

The first family here that is multivariate and **not elliptical**, and
therefore the second real test of that layer: there is no location and
scale to separate, the support is a simplex rather than a Euclidean
space, and the covariance is singular by construction because the
coordinates sum to one.

The parametrisation follows the same design as the multivariate
gaussian's. The constrained object — here a point of the simplex rather
than a positive definite matrix — is carried by a parameters7 parameter
and **flattened into scalars** with identity links, so every generic of
the package indexes it as it always did. The shapes are \\\alpha =
\phi\mu\\.

**Density:** \$\$f(y) = \dfrac{\Gamma(\phi)}{\prod_j\Gamma(\alpha_j)}
\prod_j y_j^{\alpha_j-1}\$\$

**Moments:** mean \\\mu\\ and \\\operatorname{Cov}(Y_i,Y_j) =
(\delta\_{ij}\mu_i - \mu_i\mu_j)/(\phi+1)\\, so \\\phi\\ is a precision:
the larger it is, the tighter the draws about the mean.

**The expected information is closed form**, which two identities make
possible. Differentiating \\\sum_j\mu_j = 1\\ once and twice shows that
the columns of \\A = \partial\mu/\partial\eta\\ sum to zero and so does
every second-derivative vector of the simplex; and \\\mathbb{E}\[\log
y_j\] = \psi(\alpha_j) - \psi(\phi)\\ makes \\\mathbb{E}\[g_j\] =
-\psi(\phi)\\ the same for every \\j\\. Every term carrying the data is
therefore a constant times one of those zero sums, and drops out.

**The marginals are Beta**, coordinate \\j\\ being
\\\mathrm{Beta}(\alpha_j, \phi-\alpha_j)\\, so
[`mv_marginal`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
returns an object rather than refusing — which is what makes this family
a useful test of that generic rather than another refusal. Several
coordinates together are again Dirichlet, but only after the remaining
mass is collapsed into a coordinate of its own, so that case is refused
rather than returned under a name that would mislead.

The distribution function and the quantile are refused by
[`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md),
as for every family of that class.

## See also

[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for the two-coordinate case seen on the line,
[`mvgaussian_distrib`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md),
[`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.html)

## Examples

``` r
d <- dirichlet_distrib(3)
d@params
#> [1] "mean_alr1" "mean_alr2" "phi"      

theta <- as.list(stats::setNames(c(0.3, -0.2, log(12)), d@params))
mv_location(d, theta)
#> [1] 0.4260125 0.2583897 0.3155978
round(mv_sigma(d, theta), 5)
#>          [,1]     [,2]     [,3]
#> [1,]  0.07017 -0.03159 -0.03858
#> [2,] -0.03159  0.05499 -0.02340
#> [3,] -0.03858 -0.02340  0.06198

# the marginals are Beta, so a panel of a pairs plot is a real object
mv_marginal(d, theta, which = 1)$theta
#> $mu
#> [1] 0.4260125
#> 
#> $phi
#> [1] 2.484907
#> 
```
