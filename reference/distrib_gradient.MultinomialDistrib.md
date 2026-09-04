# Multinomial Score

Computes the first derivatives of the multinomial log-mass with respect
to the simplex's free values, one value per observation, in closed form.
With \\A = \partial p/\partial\eta\\ the chart's own Jacobian,
\$\$\dfrac{\partial\ell}{\partial\eta_k} = \sum\_{j=1}^{p}
\dfrac{y_j}{p_j}A\_{jk}.\$\$ The factorial term of the mass carries no
parameter, so it contributes nothing, and the data enter only through
the counts themselves.

With `scale = "link"` the generic applies the chain rule for the links
the family carries. Every parameter rides the identity here, its free
value being unconstrained already, so the two scales coincide.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- y:

  A numeric matrix with one row per observation and \\p\\ columns, each
  row on the support. A single observation may be given as a plain
  vector.

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per free value, in the order
`distrib@params` gives, each of length `nrow(y)`.

## Notation

\\\ell\\ is the log-mass of one observation, \\p\\ the probability
vector, \\n\\ the trial count, \\\eta\\ the free vector of the simplex
chart and \\A = \partial p/\partial\eta\\ its Jacobian.

## See also

[`distrib_hessian.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MultinomialDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MultinomialDistrib.md)
for their expectation,
[`mn_parts()`](https://statmodels7.github.io/distributions7/reference/mn_parts.md)
for the shared pieces, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
set.seed(1)
Y <- distrib_rng(d, 4, th)
g <- distrib_gradient(d, Y, th)
names(g)
#> [1] "probs_alr1" "probs_alr2"

# It is the derivative of the log-mass, so numDeriv on the summed log-mass
# reproduces the summed score.
fn <- function(v)
  sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
rbind(numeric = numDeriv::grad(fn, unlist(th)),
      closed = vapply(g, sum, numeric(1)))
#>         probs_alr1 probs_alr2
#> numeric -0.5202503   1.832207
#> closed  -0.5202503   1.832207

# The score has mean zero over the support: the first Bartlett identity,
# here an exact sum.
supp <- mv_support(d, th)
w <- distrib_pdf(d, supp, th)
vapply(distrib_gradient(d, supp, th), function(v) sum(w * v), numeric(1))
#>   probs_alr1   probs_alr2 
#> 6.722053e-17 2.762547e-16 
```
