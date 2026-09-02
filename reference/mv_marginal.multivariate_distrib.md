# No Marginal Without a Closed Form

Signals an error. Integrating the other coordinates out has no general
closed form, and a numerical marginal would be a different object under
the same name: a quadrature's answer is a grid of numbers, not a
distribution another generic can be asked of. A family whose marginals
are known registers its own method, and all four that ship do.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- theta:

  A named list of parameters. Not examined: the error is raised before
  it is read.

- which:

  An integer vector of coordinates, already validated by the generic.
  Not examined.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Never returns: it always signals an error naming the family.

## See also

[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the generic and what each shipped family gives, and
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
which exists exactly where the marginals do.

## Examples

``` r
# Every shipped family registers a method, so reach the refusal directly.
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0)
base <- S7::method(mv_marginal, multivariate_distrib)
try(base(d, distributions7:::align_theta(d, theta), 1L))
#> Error : mv_marginal() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': integrating out the other coordinates has no closed form for this family, and a numerical marginal would not be the same object.

# The family's own method answers.
mv_marginal(d, theta, 1)$distrib@n_dim
#> [1] 1
```
