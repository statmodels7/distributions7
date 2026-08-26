# Beta-Binomial Fourth-Order Derivatives in Its Shapes

Computes the five distinct fourth derivatives of the beta-binomial
log-mass in the two shapes, in closed form, by the construction
[`distrib_deriv3.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom2Distrib.md)
describes carried one order further: the same sum of terms with
\\\psi_3\\, the third polygamma, in place of \\\log\Gamma\\.

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

A named list of five numeric vectors, `alpha_alpha_alpha_alpha`,
`alpha_alpha_alpha_beta`, `alpha_alpha_beta_beta`,
`alpha_beta_beta_beta` and `beta_beta_beta_beta`, each of length
`max(length(y), length(alpha), length(beta))`.

## See also

[`distrib_deriv3.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom2Distrib.md)
for the order below,
[`betabinom2_derivs()`](https://statmodels7.github.io/distributions7/reference/betabinom2_derivs.md)
for the assembly, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
th <- list(alpha = 2, beta = 3)
d4 <- distrib_deriv4(d, 0:10, th)
names(d4)
#> [1] "alpha_alpha_alpha_alpha" "alpha_alpha_alpha_beta" 
#> [3] "alpha_alpha_beta_beta"   "alpha_beta_beta_beta"   
#> [5] "beta_beta_beta_beta"    

# A central difference of the third order reproduces a mixed component.
eps <- 1e-5
up <- distrib_deriv3(d, 0:10, list(alpha = 2, beta = 3 + eps))$alpha_alpha_beta
dn <- distrib_deriv3(d, 0:10, list(alpha = 2, beta = 3 - eps))$alpha_alpha_beta
all.equal((up - dn) / (2 * eps), d4$alpha_alpha_beta_beta, tolerance = 1e-5)
#> [1] TRUE

# The expected branch is the mass-weighted sum, exactly.
w <- distrib_pdf(d, 0:10, th)
rbind(expected = vapply(distrib_deriv4(d, 0, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      summed = vapply(d4, function(v) sum(w * v), numeric(1)))
#>          alpha_alpha_alpha_alpha alpha_alpha_alpha_beta alpha_alpha_beta_beta
#> expected              -0.4148922             0.02077335            0.02077335
#> summed                -0.4148922             0.02077335            0.02077335
#>          alpha_beta_beta_beta beta_beta_beta_beta
#> expected           0.02077335          -0.0908373
#> summed             0.02077335          -0.0908373
```
