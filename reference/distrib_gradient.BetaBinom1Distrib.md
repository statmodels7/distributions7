# Beta-Binomial Score

Computes the first derivatives of the beta-binomial log-mass with
respect to the mean proportion \\\mu\\ and the dispersion \\\sigma\\,
one value per observation, in closed form. The parameters enter only
through the two beta shapes, where each derivative is a difference of
digammas: \$\$\dfrac{\partial \ell}{\partial \alpha} = \psi(y+\alpha) -
\psi(\alpha) - \psi(n+S) + \psi(S), \qquad S = \alpha + \beta,\$\$ and
likewise in \\\beta\\ with \\n - y\\ in place of \\y\\. The reported
components follow by the chain rule of \\(\alpha, \beta) = (\mu/\sigma,
(1-\mu)/\sigma)\\.

Each digamma difference cancels to leading order as the shapes grow, so
the compiled kernel forms it as a sum of reciprocals, the shifts being
integers. The location component then holds its accuracy to a
concentration of \\10^{15}\\, where the direct difference has lost every
digit.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale. The arithmetic runs in a compiled kernel decomposed
over the elements of the output, so the result does not depend on the
thread count.

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

A named list of two numeric vectors, `mu` and `sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \in (0,1)\\ the mean
proportion, \\\sigma \> 0\\ the dispersion, \\n\\ the trial count and
\\\psi\\ the digamma function,
[`digamma()`](https://rdrr.io/r/base/Special.html) in R.

## See also

[`distrib_hessian.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaBinom1Distrib.md)
for their expectation,
[`distrib_gradient.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom2Distrib.md)
for the same quantity in the shapes, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)
g <- distrib_gradient(d, 0:10, th)

# It is the derivative of the log-mass, so a central difference of it
# reproduces the mu component.
eps <- 1e-6
all.equal((distrib_pdf(d, 0:10, list(mu = 0.3 + eps, sigma = 0.5), log = TRUE) -
           distrib_pdf(d, 0:10, list(mu = 0.3 - eps, sigma = 0.5), log = TRUE)) /
            (2 * eps), g$mu, tolerance = 1e-6)
#> [1] TRUE

# The score has mean zero over the support: the first Bartlett identity.
w <- distrib_pdf(d, 0:10, th)
vapply(g, function(v) sum(w * v), numeric(1))
#>            mu         sigma 
#> -1.098080e-15 -2.463307e-16 

# At a tiny dispersion the family is a binomial, and the mu score is the
# binomial one, y/mu - (n-y)/(1-mu), evaluated here at y = 3.
c(kernel = distrib_gradient(d, 3, list(mu = 0.4, sigma = 1e-12))$mu,
  binomial = 3 / 0.4 - 7 / 0.6)
#>    kernel  binomial 
#> -4.166667 -4.166667 
```
