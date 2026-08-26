# Beta-Binomial Fourth-Order Derivatives in Mean and Dispersion

Computes the five distinct fourth derivatives of the beta-binomial
log-mass in \\\mu\\ and \\\sigma\\, **in closed form**, by the
construction
[`distrib_deriv3.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom1Distrib.md)
describes carried one order further: the shape parametrization's fourth
derivatives, each a difference of \\\psi^{(3)}\\, carried through the
map \\(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)\\ by the partition
sum of
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md).

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead. That is the one place on this page where `approx` and `nsim`
are read.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

- y:

  A numeric vector of counts in \\\\0, \dots, n\\\\. With
  `expected = TRUE` only its length is read.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `mu` must lie in \\(0, 1)\\ and `sigma` be strictly
  positive.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the value at the data. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  One of `"integrate"` (the default here), `"bartlett"`, `"mc"` or
  `"opg"`, the strategy
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  uses. Read only when `expected = TRUE`.

- nsim:

  A single positive integer, the sample size when `approx = "mc"`. Read
  only when `expected = TRUE`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma`,
`mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
`sigma_sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \in (0,1)\\ the mean
proportion, \\\sigma \> 0\\ the dispersion and \\n\\ the trial count.
\\\alpha\\ and \\\beta\\ are the two beta shapes the family is written
in internally.

## See also

[`distrib_deriv3.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom1Distrib.md)
for the order below,
[`distrib_hessian.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom1Distrib.md)
for the second order,
[`betabinom1_components()`](https://statmodels7.github.io/distributions7/reference/betabinom1_components.md)
for the assembly, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)
d4 <- distrib_deriv4(d, 0:10, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# A central difference of the third order reproduces a mixed component.
eps <- 1e-5
up <- distrib_deriv3(d, 0:10, list(mu = 0.3, sigma = 0.5 + eps))$mu_mu_sigma
dn <- distrib_deriv3(d, 0:10, list(mu = 0.3, sigma = 0.5 - eps))$mu_mu_sigma
all.equal((up - dn) / (2 * eps), d4$mu_mu_sigma_sigma, tolerance = 1e-5)
#> [1] TRUE

# The expected branch is a numerical expectation, and the mass-weighted
# sum over the support reaches it.
w <- distrib_pdf(d, 0:10, th)
rbind(expected = vapply(distrib_deriv4(d, 0, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      summed = vapply(d4, function(v) sum(w * v), numeric(1)))
#>          mu_mu_mu_mu mu_mu_mu_sigma mu_mu_sigma_sigma mu_sigma_sigma_sigma
#> expected   -582.8087      -6.091053         -27.80562              4.46073
#> summed     -582.8087      -6.091053         -27.80562              4.46073
#>          sigma_sigma_sigma_sigma
#> expected                -47.0173
#> summed                  -47.0173
```
