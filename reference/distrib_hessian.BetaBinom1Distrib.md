# Beta-Binomial Observed Hessian

Computes the three distinct second derivatives of the beta-binomial
log-mass with respect to \\\mu\\ and \\\sigma\\, one value per
observation, in closed form. In the shapes each second derivative is a
difference of trigammas, \$\$\dfrac{\partial^2 \ell}{\partial \alpha^2}
= \psi_1(y+\alpha) - \psi_1(\alpha) - \psi_1(n+S) + \psi_1(S), \qquad S
= \alpha + \beta,\$\$ and the **mixed** shape component carries only the
\\S\\ part, \\-\psi_1(n+S) + \psi_1(S)\\, the two shapes entering the
mass separately otherwise. The reported components follow by the
two-variable chain rule of \\(\alpha, \beta) = (\mu/\sigma,
(1-\mu)/\sigma)\\, whose own second derivatives contribute the rest.

The arithmetic runs in a compiled kernel decomposed over the elements of
the output, so the result does not depend on the thread count.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

- y:

  A numeric vector of counts in \\\\0, \dots, n\\\\.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `mu` must lie in \\(0, 1)\\ and `sigma` be strictly
  positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_sigma` and
`sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`. The three name the distinct
entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \in (0,1)\\ the mean
proportion, \\\sigma \> 0\\ the dispersion, \\n\\ the trial count and
\\\psi_1\\ the trigamma function,
[`trigamma()`](https://rdrr.io/r/base/Special.html) in R.

## See also

[`distrib_gradient.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom1Distrib.md)
for the score,
[`distrib_expected_hessian.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaBinom1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)
h <- distrib_hessian(d, 0:10, th)
names(h)
#> [1] "mu_mu"       "mu_sigma"    "sigma_sigma"

# A central difference of the score reproduces the pure-mu component.
eps <- 1e-5
up <- distrib_gradient(d, 0:10, list(mu = 0.3 + eps, sigma = 0.5))$mu
dn <- distrib_gradient(d, 0:10, list(mu = 0.3 - eps, sigma = 0.5))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#> [1] TRUE

# The mass-weighted sum over the support is the expected Hessian.
w <- distrib_pdf(d, 0:10, th)
rbind(summed = vapply(h, function(v) sum(w * v), numeric(1)),
      expected = vapply(distrib_expected_hessian(d, 0:10, th),
                        function(v) v[1], numeric(1)))
#>              mu_mu  mu_sigma sigma_sigma
#> summed   -12.99319 0.9716737   -1.208422
#> expected -12.99319 0.9716737   -1.208422
```
