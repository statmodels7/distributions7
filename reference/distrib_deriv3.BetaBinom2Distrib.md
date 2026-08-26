# Beta-Binomial Third-Order Derivatives in Its Shapes

Computes the four distinct third derivatives of the beta-binomial
log-mass in the two shapes, in closed form. The log-mass is a sum of
log-gamma terms, so the third derivative is the same sum with
\\\psi_2\\, the second polygamma, in place of \\\log\Gamma\\. A
component naming both shapes keeps only the terms in \\\alpha+\beta\\
and is free of the data.

With `expected = TRUE` the expectation is an **exact finite sum** over
the support, so `approx` and `nsim` are ignored on both branches.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

- y:

  A numeric vector of counts in \\\\0, \dots, n\\\\. With
  `expected = TRUE` only its length is read.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the exact expectation under the model
  is returned in place of the value at the data. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored here, both branches being exact. Accepted so that the
  signature matches the generic's.

- nsim:

  Ignored here, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of four numeric vectors, `alpha_alpha_alpha`,
`alpha_alpha_beta`, `alpha_beta_beta` and `beta_beta_beta`, each of
length `max(length(y), length(alpha), length(beta))`.

## See also

[`distrib_hessian.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom2Distrib.md)
for the order below,
[`distrib_deriv4.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaBinom2Distrib.md)
for the order above,
[`betabinom2_derivs()`](https://statmodels7.github.io/distributions7/reference/betabinom2_derivs.md)
for the assembly, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
th <- list(alpha = 2, beta = 3)
d3 <- distrib_deriv3(d, 0:10, th)
names(d3)
#> [1] "alpha_alpha_alpha" "alpha_alpha_beta"  "alpha_beta_beta"  
#> [4] "beta_beta_beta"   

# A central difference of the Hessian reproduces the pure-alpha component.
eps <- 1e-5
up <- distrib_hessian(d, 0:10, list(alpha = 2 + eps, beta = 3))$alpha_alpha
dn <- distrib_hessian(d, 0:10, list(alpha = 2 - eps, beta = 3))$alpha_alpha
all.equal((up - dn) / (2 * eps), d3$alpha_alpha_alpha, tolerance = 1e-6)
#> [1] TRUE

# The expected branch is the mass-weighted sum, exactly.
w <- distrib_pdf(d, 0:10, th)
rbind(expected = vapply(distrib_deriv3(d, 0, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      summed = vapply(d3, function(v) sum(w * v), numeric(1)))
#>          alpha_alpha_alpha alpha_alpha_beta alpha_beta_beta beta_beta_beta
#> expected         0.2872239      -0.04403913     -0.04403913     0.09004371
#> summed           0.2872239      -0.04403913     -0.04403913     0.09004371
```
