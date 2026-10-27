# Beta Expected Hessian in the Shapes

Returns the same three numbers as
[`distrib_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md),
because the observed Hessian is free of the response and so equals its
own expectation: \$\$\mathbb{E}\left\[\ell^{(\alpha\alpha)}\right\] =
\psi'(\alpha+\beta) - \psi'(\alpha), \qquad
\mathbb{E}\left\[\ell^{(\alpha\beta)}\right\] = \psi'(\alpha+\beta),
\qquad \mathbb{E}\left\[\ell^{(\beta\beta)}\right\] =
\psi'(\alpha+\beta) - \psi'(\beta).\$\$ The expectation is exact, so the
method does not average, integrate or simulate, and it reads `y` only
for its length.

The consequence for fitting is that Fisher scoring and Newton's method
take the same step on the parameter scale, because both invert the same
matrix. On the link scale they differ, because the chain rule there adds
a term in the score, which does read the data.

The mixed entry \\\psi'(\alpha+\beta)\\ is positive and never zero, so
the two shapes are not orthogonal at any parameter setting.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  a fixed number of observations, set per kernel, the kernel runs
  sequentially. Defaults to `1L`.

## Value

A named list of three numeric vectors, `alpha_alpha`, `beta_beta` and
`alpha_beta`, in that order, each of length
`max(length(y), length(alpha), length(beta))`.

## Notation

The expected information is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the observed information under the model. The beta is
a regular family, so the second Bartlett identity holds and this equals
the variance of the score.

## See also

[`distrib_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md),
which returns the same numbers;
[`distrib_expected_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta1Distrib.md),
where the same law in the mean and the precision does have an observed
Hessian that reads the data;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md);
and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)
lapply(distrib_expected_hessian(d, y, th), unique)
#> $alpha_alpha
#> [1] -0.4913889
#> 
#> $beta_beta
#> [1] -0.06777778
#> 
#> $alpha_beta
#> [1] 0.1535452
#> 

# Identical to the observed Hessian, so the two fitting methods coincide.
identical(distrib_expected_hessian(d, y, th), distrib_hessian(d, y, th))
#> [1] TRUE

set.seed(8)
z <- distrib_rng(d, 2000, th)
rbind(newton = coef(fit_distrib(d, z, method = optimizers7::newton())),
      fisher = coef(fit_distrib(d, z, method = fisher_scoring())))
#>           alpha     beta
#> newton 1.927355 4.827894
#> fisher 1.927355 4.827894

# The mixed entry is psi'(alpha + beta), positive at every setting, so the
# two shapes are never orthogonal.
distrib_expected_hessian(d, 0.5, th)$alpha_beta
#> [1] 0.1535452
```
