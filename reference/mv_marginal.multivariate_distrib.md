# Marginal of a Family Without a Closed Form

Signals an error. Integrating the other coordinates out has no general
closed form, and a numerical marginal would be a grid of values and not
a distribution object to which the generics apply. A family whose
marginals are known registers its own method, as every shipped family
does.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- theta:

  A named list of parameters. Not read: the error is signalled first.

- which:

  An integer vector of coordinates, already validated by the generic.
  Not read.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

Does not return: it signals an error naming the family.

## See also

[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the generic and what each shipped family returns, and
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
whose panels are marginals.

## Examples

``` r
# Every shipped family registers a method, so the method of the class is
# called directly.
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0)
base <- S7::method(mv_marginal, multivariate_distrib)
try(base(d, distributions7:::align_theta(d, theta), 1L))
#> Error : mv_marginal() is not defined for 'multivariate gaussian [2d, sigma=log_cholesky]': integrating out the other coordinates has no closed form for this family, and a numerical marginal would not be the same object.

# The family's own method returns the marginal.
mv_marginal(d, theta, 1)$distrib@n_dim
#> [1] 1
```
