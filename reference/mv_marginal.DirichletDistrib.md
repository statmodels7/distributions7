# Dirichlet Marginal

Returns the law of one coordinate, which is \\\mathrm{Beta}(\alpha_j,
\phi - \alpha_j)\\. In this package's mean and precision parametrization
of the beta that is mean \\\mu_j\\ and precision \\\phi\\, so every
marginal carries the concentration of the joint law, as every marginal
of the multivariate t carries its degrees of freedom.

The marginal of one coordinate is a univariate beta distribution, so,
unlike the marginals of the elliptical families, it has a distribution
function, a quantile and the moments of a univariate family. A vector of
several coordinates signals an error: a sub-vector of a Dirichlet is
again Dirichlet only once the remaining mass is collapsed into a
coordinate of its own, which is a different object from the sub-vector.
For the same reason
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
whose panels are marginals of two coordinates, signals that error for
this family.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

- which:

  A single integer in \\1, \dots, p\\, the coordinate wanted. A vector
  of two or more distinct coordinates signals an error explaining why;
  an empty, duplicated, non-integer or out-of-range `which` is rejected
  earlier by
  [`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `distrib`, a
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
object, and `theta`, a list with the marginal's `mu` and `phi`.

## See also

[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for the family returned,
[`mv_sigma.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.DirichletDistrib.md)
for the coordinate variances, and
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
m <- mv_marginal(d, th, which = 1)
m$distrib
#> Distribution: Beta1
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu  (mean)               | Link: logit      | Domain: (0, 1)
#>   phi (precision)          | Link: log        | Domain: (0, Inf)
m$theta
#> $mu
#> [1] 0.4260125
#> 
#> $phi
#> [1] 12
#> 

# The marginal's mean and variance are the first coordinate's.
rbind(marginal = c(mean(m$distrib, m$theta), variance(m$distrib, m$theta)),
      joint = c(mv_location(d, th)[1], mv_sigma(d, th)[1, 1]))
#>               [,1]       [,2]
#> marginal 0.4260125 0.01880968
#> joint    0.4260125 0.01880968

# Every marginal carries the same concentration.
vapply(1:3, function(j) mv_marginal(d, th, which = j)$theta$phi, numeric(1))
#> [1] 12 12 12
```
