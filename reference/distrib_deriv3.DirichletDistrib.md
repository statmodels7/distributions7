# Dirichlet Third Derivatives

Computes every third derivative of the log-density in the parameters, in
closed form. The log-density is \$\$\ell = \log\Gamma(\phi) - \sum_j
\log\Gamma(\alpha_j) + \sum_j(\alpha_j - 1)\log y_j, \qquad \alpha =
\phi\\\mu(\eta),\$\$ so every term depends on one coordinate of the
simplex and the chain rule reduces to a univariate partition sum per
coordinate, run by
[`chain_univariate()`](https://statmodels7.github.io/distributions7/reference/chain_univariate.md)
over the arrays that
[`dirichlet_map_tensors()`](https://statmodels7.github.io/distributions7/reference/dirichlet_map_tensors.md)
supplies. The concentration also enters directly, through
\\\log\Gamma(\phi)\\, which contributes only to the component whose
indices all name it.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- y:

  A numeric matrix with one row per observation and one column per
  coordinate, where each row is a point of the open simplex. A row with
  a zero coordinate is outside the support and gives a non-finite
  component through \\\log y_j\\.

- theta:

  A named list of parameters, where each component is a single number:
  `mean_alr1`, ..., `mean_alr(p-1)` and `phi`.

- expected:

  Logical of length 1. When `TRUE` the expectation is returned, in
  closed form (see Details). Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`. The method returns the
  parameter scale and the generic applies the link scale to its result.

- approx:

  Not read: the expectation is in closed form. It is present so that the
  signature matches the generic's.

- nsim:

  Not read. It is present so that the signature matches the generic's.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 3)`. At \\p = 3\\ there are ten components.
With `expected = TRUE` every vector is constant.

## Details

The shape vector is bilinear in \\\phi\\ and \\\mu\\, so the arrays of
the map are zero for a component that names \\\phi\\ twice or more.

On the additive log-ratio chart of
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
the term \\-\log\Gamma(\alpha_j)\\ is split as \\\log\alpha_j -
\log\Gamma(1 + \alpha_j)\\. The first part equals \\\log\phi + \eta_j -
\mathrm{LSE}(\eta)\\, with \\\eta_p = 0\\, and its derivatives are added
in closed form: \\p(-1)^{m-1}(m-1)!/\phi^m\\ on the component of order
\\m\\ in \\\phi\\ alone, and \\-p\\ times a cumulant of the categorical
law (see
[`mn_cumulant()`](https://statmodels7.github.io/distributions7/reference/mn_cumulant.md))
on a component that names only the free values of the mean. The second
part is regular at a small shape, so the partition sum runs over
\\\psi^{(m-1)}(1 + \alpha_j)\\. Applied to \\\log\Gamma(\alpha_j)\\
directly, the chain rule would cancel at a small shape: terms of order
one would cancel and leave a result of the order of \\\alpha_j\\. On
another chart a coordinate with \\\alpha_j \< 1\\ is scaled by
\\\alpha_j\\ through \\\alpha^m\psi^{(m-1)}(\alpha) = (-1)^m(m-1)! +
\alpha^m\psi^{(m-1)}(1 + \alpha)\\, so that \\\psi^{(m-1)}(\alpha_j)\\
does not overflow.

With `expected = TRUE` the expectation is exact. The response enters
only through \\\sum_j(\alpha_j - 1)\log y_j\\, linearly, and
\\\mathrm{E}\log y_j = \psi(\alpha_j) - \psi(\phi)\\.

## Notation

\\\alpha\\ is the shape vector, \\\phi\\ the concentration, \\\mu\\ the
mean on the simplex, \\\eta\\ its free vector, \\p\\ the number of
coordinates, \\\Gamma\\ the gamma function, \\\psi\\ the digamma
function and \\\ell\\ the log-density of one observation.

## See also

[`distrib_deriv4.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.DirichletDistrib.md)
for the next order,
[`distrib_hessian.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.DirichletDistrib.md)
for the second,
[`dirichlet_higher()`](https://statmodels7.github.io/distributions7/reference/dirichlet_higher.md)
for the shared engine, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
set.seed(1)
y <- distrib_rng(d, 5, theta)

d3 <- distrib_deriv3(d, y, theta)
names(d3)
#>  [1] "mean_alr1_mean_alr1_mean_alr1" "mean_alr1_mean_alr1_mean_alr2"
#>  [3] "mean_alr1_mean_alr1_phi"       "mean_alr1_mean_alr2_mean_alr2"
#>  [5] "mean_alr1_mean_alr2_phi"       "mean_alr1_phi_phi"            
#>  [7] "mean_alr2_mean_alr2_mean_alr2" "mean_alr2_mean_alr2_phi"      
#>  [9] "mean_alr2_phi_phi"             "phi_phi_phi"                  

# The result is compared with a central difference of the Hessian in phi.
h <- 1e-4
tp <- theta; tp$phi <- tp$phi + h
tm <- theta; tm$phi <- tm$phi - h
c(exact = sum(d3[["mean_alr1_mean_alr1_phi"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["mean_alr1_mean_alr1"]]) -
             sum(distrib_hessian(d, y, tm)[["mean_alr1_mean_alr1"]])) / (2 * h))
#>     exact   stencil 
#> -1.096431 -1.096431 

# The pure-phi component carries no response: it is psi''(phi) minus the
# mu_j^3-weighted psi''(phi mu_j), the same at every observation.
mu <- as.numeric(parameters7::param_value(d@param, c(0.3, -0.2)))
c(component = d3[["phi_phi_phi"]][1],
  direct = psigamma(8, 2) - sum(mu^3 * psigamma(8 * mu, 2)))
#>   component      direct 
#> 0.004874712 0.004874712 

# The expectation is exact, so it does not depend on the seed.
distrib_deriv3(d, y, theta, expected = TRUE)[["phi_phi_phi"]][1]
#> [1] 0.004874712
```
