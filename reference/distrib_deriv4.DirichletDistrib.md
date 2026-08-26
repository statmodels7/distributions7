# Dirichlet Fourth Derivatives

Computes every fourth derivative of the log-density in the parameters,
in closed form, by the same construction as
[`distrib_deriv3.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.DirichletDistrib.md)
one order up: a univariate partition sum per coordinate of the simplex,
run by
[`chain_univariate()`](https://statmodels7.github.io/distributions7/reference/chain_univariate.md)
over the arrays
[`dirichlet_map_tensors()`](https://statmodels7.github.io/distributions7/reference/dirichlet_map_tensors.md)
supplies, plus the direct \\\log\Gamma(\phi)\\ term on the component all
of whose indices name the concentration.

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
`deriv_names(distrib@params, 4)`. At \\p = 3\\ there are fifteen
components. With `expected = TRUE` every vector is constant.

## Details

The license for this order is that the SAME assembly run at orders one
and two reproduces the hand-written score and information, which are
derived separately and are already under
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md).
Agreement at the orders that can be checked is what authorizes the order
that cannot.

With `expected = TRUE` the expectation is taken by sampling and carries
Monte Carlo error of order `nsim^(-1/2)`.

## Notation

\\\alpha\\ is the shape vector, \\\phi\\ the concentration, \\\mu\\ the
mean on the simplex, \\\eta\\ its free vector and \\\ell\\ the
log-density of one observation.

## See also

[`distrib_deriv3.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.DirichletDistrib.md)
for the order below,
[`dirichlet_higher()`](https://statmodels7.github.io/distributions7/reference/dirichlet_higher.md)
for the shared engine, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
set.seed(1)
y <- distrib_rng(d, 5, theta)

d4 <- distrib_deriv4(d, y, theta)
length(d4)
#> [1] 15

# Against one stencil on the analytic third order.
h <- 1e-4
tp <- theta; tp$phi <- tp$phi + h
tm <- theta; tm$phi <- tm$phi - h
c(exact = sum(d4[["mean_alr1_mean_alr1_phi_phi"]]),
  stencil = (sum(distrib_deriv3(d, y, tp)[["mean_alr1_mean_alr1_phi"]]) -
             sum(distrib_deriv3(d, y, tm)[["mean_alr1_mean_alr1_phi"]])) /
            (2 * h))
#>        exact      stencil 
#> -0.006042767 -0.006042767 
```
