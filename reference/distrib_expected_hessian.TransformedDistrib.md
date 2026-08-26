# Transformed Expected Information

Returns the parent's expected Hessian, EXACTLY. Since \\\ell_Y(\theta;
y) = \ell_X(\theta; g^{-1}(y)) + \log\lvert J(y)\rvert\\ and the
Jacobian does not depend on \\\theta\\, the expectation of the
transformed model's second derivative is the parent's expectation
reparametrized by the change of variables, which is the same number. No
approximation runs and no Monte Carlo is needed, whatever `approx` says.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- y:

  A numeric vector of observations on the transformed scale. Read only
  through the parent, which uses its length.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- approx:

  Ignored: the delegation is exact. Present so that the signature
  matches the generic's. Note that it is not forwarded, so a parent
  whose OWN expected Hessian is approximate takes its own default.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

This is the strongest of the three delegations: the score and the
observed Hessian are the parent's AT A POINT, and this one is the
parent's as a whole. A family that is a transformation of a family with
a closed-form information therefore has one too, at no cost.

## Notation

\\g\\ is the transformation, \\J\\ the Jacobian of its inverse and
\\\ell\\ a log-density.

## See also

[`distrib_hessian.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TransformedDistrib.md)
for the observed matrix,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
whose Fisher scoring inverts this, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)
set.seed(2)
y <- distrib_rng(d, 40, theta)

EH <- distrib_expected_hessian(d, y, theta)
vapply(EH, function(z) z[1], numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>     -1.5625     -3.1250      0.0000 

# It is the parent's, at the preimage.
identical(EH, distrib_expected_hessian(gaussian1_distrib(), log(y), theta))
#> [1] TRUE

# So a lognormal built this way inherits the gaussian's closed form, whose
# mean block is -1 / sigma^2.
c(reported = EH$mu_mu[1], theory = -1 / 0.8^2)
#> reported   theory 
#>  -1.5625  -1.5625 
```
