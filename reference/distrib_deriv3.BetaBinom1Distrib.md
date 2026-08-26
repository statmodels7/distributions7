# Beta-Binomial Third-Order Derivatives in Mean and Dispersion

Computes the four distinct third derivatives of the beta-binomial
log-mass in \\\mu\\ and \\\sigma\\, **in closed form**. The shape
parametrization carries closed derivatives at every order, each a
difference of polygammas, and this parametrization is that one at
\\\alpha = \mu/\sigma\\ and \\\beta = (1-\mu)/\sigma\\, so the Faa di
Bruno partition sum of
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
over the map delivers them. Both shapes are linear in \\\mu\\ at fixed
\\\sigma\\, so every partial of the map carrying two or more \\\mu\\
vanishes and the sum is short.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead. That is the one place on this page where `approx` and `nsim`
are read; on the observed branch both are ignored.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`. The four name the distinct
entries of a symmetric third-order array over two parameters.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \in (0,1)\\ the mean
proportion, \\\sigma \> 0\\ the dispersion and \\n\\ the trial count.
\\\alpha\\ and \\\beta\\ are the two beta shapes the family is written
in internally.

## See also

[`distrib_hessian.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom1Distrib.md)
for the order below,
[`distrib_deriv4.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaBinom1Distrib.md)
for the order above,
[`betabinom1_components()`](https://statmodels7.github.io/distributions7/reference/betabinom1_components.md)
for the assembly, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)
d3 <- distrib_deriv3(d, 0:10, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# A central difference of the Hessian reproduces the pure-mu component,
# which is what says the partition sum over the map is right.
eps <- 1e-5
up <- distrib_hessian(d, 0:10, list(mu = 0.3 + eps, sigma = 0.5))$mu_mu
dn <- distrib_hessian(d, 0:10, list(mu = 0.3 - eps, sigma = 0.5))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#> [1] TRUE

# The expected branch is a numerical expectation of the observed one, and
# the mass-weighted sum over the support reaches it.
w <- distrib_pdf(d, 0:10, th)
rbind(expected = vapply(distrib_deriv3(d, 0, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      summed = vapply(d3, function(v) sum(w * v), numeric(1)))
#>          mu_mu_mu mu_mu_sigma mu_sigma_sigma sigma_sigma_sigma
#> expected  49.9161    7.067264      -1.768251          6.964833
#> summed    49.9161    7.067264      -1.768251          6.964833
```
