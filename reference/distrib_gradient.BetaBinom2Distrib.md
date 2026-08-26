# Beta-Binomial Score in Its Shapes

Computes the first derivatives of the beta-binomial log-mass with
respect to the two shapes, one value per observation, in closed form:
\$\$\dfrac{\partial\ell}{\partial\alpha} = \psi(y+\alpha) -
\psi(\alpha) - \psi(n+\alpha+\beta) + \psi(\alpha+\beta),\$\$ and the
same with \\n-y\\ and \\\beta\\ in place of \\y\\ and \\\alpha\\. The
two share the term in \\\alpha+\beta\\, which is the only part a mixed
second derivative keeps.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale; the transformation happens in the generic.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

- y:

  A numeric vector of counts in \\\\0, \dots, n\\\\.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `alpha` and `beta`, each of length
`max(length(y), length(alpha), length(beta))`.

## Notation

\\\ell\\ is the log-mass of one observation, \\\alpha, \beta \> 0\\ the
two beta shapes, \\n\\ the trial count and \\\psi\\ the digamma
function, [`digamma()`](https://rdrr.io/r/base/Special.html) in R.

## Large concentration

Each digamma difference cancels to leading order as the shapes grow at a
fixed ratio, and this method writes it out directly. Measured at \\n =
10\\, \\y = 3\\ and a mean proportion of 0.4, the relative error against
the limiting \\-5/(2S)\\ with \\S = \alpha+\beta\\ is \\6\times10^{-8}\\
at \\S = 10^8\\, \\4\times10^{-4}\\ at \\10^{12}\\ and 1.8 at
\\10^{15}\\. Where a concentration of that size is reachable, use
[`distrib_gradient.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom1Distrib.md),
whose compiled kernel forms the same difference as a sum of reciprocals
and holds nine figures there.

## See also

[`distrib_hessian.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom2Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BetaBinom2Distrib.md)
for their expectation,
[`distrib_gradient.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom1Distrib.md)
for the same quantity in the mean and dispersion, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
th <- list(alpha = 2, beta = 3)
g <- distrib_gradient(d, 0:10, th)

# The alpha component, written out.
ds <- digamma(5) - digamma(15)
all.equal(g$alpha, digamma(0:10 + 2) - digamma(2) + ds)
#> [1] TRUE

# It is the derivative of the log-mass, so a central difference reproduces
# it.
eps <- 1e-6
all.equal((distrib_pdf(d, 0:10, list(alpha = 2 + eps, beta = 3), log = TRUE) -
           distrib_pdf(d, 0:10, list(alpha = 2 - eps, beta = 3), log = TRUE)) /
            (2 * eps), g$alpha, tolerance = 1e-6)
#> [1] TRUE

# The score has mean zero over the support: the first Bartlett identity.
w <- distrib_pdf(d, 0:10, th)
vapply(g, function(v) sum(w * v), numeric(1))
#>        alpha         beta 
#> 1.491862e-16 2.107689e-16 
```
