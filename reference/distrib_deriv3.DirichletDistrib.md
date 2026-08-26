# Dirichlet Third Derivatives

Computes every third derivative of the log-density in the parameters, in
closed form. The log-density is \$\$\ell = \log\Gamma(\phi) - \sum_j
\log\Gamma(\alpha_j) + \sum_j(\alpha_j - 1)\log y_j, \qquad \alpha =
\phi\\\mu(\eta),\$\$ so every term depends on ONE coordinate of the
simplex and the chain rule collapses to a univariate partition sum per
coordinate, run by
[`chain_univariate()`](https://statmodels7.github.io/distributions7/reference/chain_univariate.md)
over the arrays
[`dirichlet_map_tensors()`](https://statmodels7.github.io/distributions7/reference/dirichlet_map_tensors.md)
supplies. The concentration also enters directly, through
\\\log\Gamma(\phi)\\, which contributes only to the component all of
whose indices name it.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- y:

  A numeric matrix with one row per observation and one column per
  coordinate, each row a point of the open simplex. A row with a zero
  coordinate is outside the support and gives a non-finite component
  through \\\log y_j\\.

- theta:

  A named list of parameters, each component a single number:
  `mean_alr1`, ..., `mean_alr(p-1)` and `phi`.

- expected:

  Logical of length 1. When `TRUE` the expectation is returned, by
  SAMPLING through
  [`mv_expected_higher()`](https://statmodels7.github.io/distributions7/reference/mv_expected_higher.md).
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- approx:

  Ignored: sampling is the only multivariate route to an expectation.
  Present so that the signature matches the generic's.

- nsim:

  The number of draws used when `expected = TRUE`. Defaults to `10000`.
  Ignored otherwise.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 3)`. At \\p = 3\\ there are ten components.
With `expected = TRUE` every vector is constant.

## Details

The shape vector is bilinear in \\\phi\\ and \\\mu\\, so a component
naming \\\phi\\ twice or more gets nothing from the map, and only the
direct \\\log\Gamma(\phi)\\ term survives there. That is what keeps the
order-3 and order-4 assemblies as cheap as the order-2 one.

With `expected = TRUE` the expectation is taken by sampling and carries
Monte Carlo error of order `nsim^(-1/2)`.

## Notation

\\\alpha\\ is the shape vector, \\\phi\\ the concentration, \\\mu\\ the
mean on the simplex, \\\eta\\ its free vector, \\\Gamma\\ the gamma
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

# Against one stencil on the analytic Hessian, which shares no algebra with
# the partition sum.
h <- 1e-4
tp <- theta; tp$phi <- tp$phi + h
tm <- theta; tm$phi <- tm$phi - h
c(exact = sum(d3[["mean_alr1_mean_alr1_phi"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["mean_alr1_mean_alr1"]]) -
             sum(distrib_hessian(d, y, tm)[["mean_alr1_mean_alr1"]])) / (2 * h))
#>     exact   stencil 
#> -1.096431 -1.096431 

# The pure-phi component gets nothing from the map beyond the shapes
# themselves: it is psi''(phi) minus the mu_j^3-weighted psi''(phi mu_j),
# and so is the same at every observation.
mu <- as.numeric(parameters7::param_value(d@param, c(0.3, -0.2)))
c(component = d3[["phi_phi_phi"]][1],
  direct = psigamma(8, 2) - sum(mu^3 * psigamma(8 * mu, 2)))
#>   component      direct 
#> 0.004874712 0.004874712 
```
