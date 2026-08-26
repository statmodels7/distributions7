# Transformed Observed Hessian

Returns the parent's observed Hessian evaluated at \\x = g^{-1}(y)\\,
for the reason the score has: the Jacobian is a constant in \\\theta\\
and differentiates away twice over.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- y:

  A numeric vector of observations on the transformed scale.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

The transformation does not depend on \\\theta\\, so \\\ell_Y(\theta; y)
= \ell_X(\theta; g^{-1}(y)) + \log\lvert J(y)\rvert\\ has a second term
that is a CONSTANT in the parameters and differentiates away. Every
derivative in \\\theta\\ is therefore the parent's own, evaluated at \\x
= g^{-1}(y)\\, and nothing is recomputed. The response derivatives are
the exception and come from the base class, the Jacobian depending on
\\y\\.

## Notation

\\g\\ is the transformation, \\J\\ the Jacobian of its inverse and
\\\ell\\ a log-density.

## See also

[`distrib_gradient.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TransformedDistrib.md)
for the first order,
[`distrib_expected_hessian.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TransformedDistrib.md)
for the expectation, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)
set.seed(2)
y <- distrib_rng(d, 40, theta)

H <- distrib_hessian(d, y, theta)
vapply(H, sum, numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>    -62.5000   -164.9888    -12.1987 

# The parent's, at the preimage.
identical(H, distrib_hessian(gaussian1_distrib(), log(y), theta))
#> [1] TRUE

# Against a numerical Hessian of the transformed log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
Hn <- numDeriv::hessian(ll, unlist(theta))
ref <- vapply(distributions7:::hess_pairs(d@params),
              function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
              numeric(1))
max(abs(vapply(H, sum, numeric(1)) - ref))
#> [1] 1.826976e-09
```
