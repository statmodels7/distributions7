# Dirichlet Marginal

Returns the law of one coordinate, which is \\\mathrm{Beta}(\alpha_j,
\phi - \alpha_j)\\. In this package's mean and precision parametrization
of the beta that is mean \\\mu_j\\ and precision \\\phi\\, so no
reparametrization is needed at all and **the concentration is shared by
every marginal**, as the multivariate t's degrees of freedom are.

This is one of the few families for which the generic returns an object
rather than signaling an error. A panel of a pairs plot is therefore a
real distribution here.

Several coordinates at once are refused. A sub-vector of a Dirichlet is
again Dirichlet, but only after the remaining mass is collapsed into a
coordinate of its own, and returning that object under this name would
mislead. The error says so.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

- which:

  A single integer in \\1, \dots, p\\, the coordinate wanted. A vector
  of length other than 1 signals an error explaining why.

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
