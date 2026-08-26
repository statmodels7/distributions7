# Beta Expected Hessian in the Shapes

Returns the same three numbers as
[`distrib_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md),
the observed Hessian being free of the response and so equal to its own
expectation: \$\$\mathbb{E}\left\[\ell^{(\alpha\alpha)}\right\] =
\psi'(\alpha+\beta) - \psi'(\alpha), \qquad
\mathbb{E}\left\[\ell^{(\alpha\beta)}\right\] = \psi'(\alpha+\beta),
\qquad \mathbb{E}\left\[\ell^{(\beta\beta)}\right\] =
\psi'(\alpha+\beta) - \psi'(\beta).\$\$ Nothing is averaged, integrated
or simulated, so `approx` and `nsim` are ignored and `y` is read only
for its length.

The consequence for fitting is that Fisher scoring and Newton's method
take the same step on the parameter scale, both inverting the same
matrix. On the **link** scale they differ, the chain rule there adding a
term in the score, which does read the data.

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
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, the expectation being exact. Accepted so that the signature
  matches the generic's, where it selects between the Bartlett,
  quadrature, Monte Carlo and outer-product routes.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `alpha_alpha`, `beta_beta` and
`alpha_beta`, in that order, each of length `length(y)`.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
beta is a regular family, so the second Bartlett identity holds and this
equals the variance of the score.

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
rbind(newton = coef(fit_distrib(d, z, method = "newton")),
      fisher = coef(fit_distrib(d, z, method = "fisher")))
#>           alpha     beta
#> newton 1.927355 4.827894
#> fisher 1.927355 4.827894

# The mixed entry is psi'(alpha + beta), positive at every setting, so the
# two shapes are never orthogonal.
distrib_expected_hessian(d, 0.5, th)$alpha_beta
#> [1] 0.1535452
```
