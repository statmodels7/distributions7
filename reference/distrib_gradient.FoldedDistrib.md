# Folded Score

Computes the first derivatives of the folded log-density,
\$\$\frac{\partial \ell}{\partial \theta_i} = w\\ s_i(x) + (1-w)\\
s_i(-x), \qquad w = \frac{f(x)}{L(x)},\$\$ the score of a TWO-COMPONENT
MIXTURE whose components are the two preimages. It needs nothing of the
fold beyond the weight: \\s\\ is the parent's own score, evaluated at
\\+x\\ and at \\-x\\.

## Arguments

- distrib:

  A `FoldedDistrib` object, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- y:

  A numeric vector of observations, non-negative.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The links are the parent's, so the two scales differ
  exactly where the parent's do.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per parameter, in `distrib@params`
order, which is the parent's order.

## Details

The weight is where the folded score departs from the parent's. Near
zero \\w\\ is one half and the two preimages contribute equally; far
into the tail of a parent centered above zero \\w\\ approaches one and
the folded score approaches the parent's.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one, \\s\\ the parent's
score, \\w\\ the weight of the positive preimage and \\\ell\\ the folded
log-density.

## See also

[`distrib_hessian.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.FoldedDistrib.md)
for the second order,
[`fold_parts()`](https://statmodels7.github.io/distributions7/reference/fold_parts.md)
for the weight, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)
set.seed(2)
y <- distrib_rng(d, 30, theta)

g <- distrib_gradient(d, y, theta)
vapply(g, sum, numeric(1))
#>        mu     sigma 
#>  4.649657 10.488692 

# Against a numerical derivative of the folded log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 2.725598e-10

# It is the mixture score, written out from the parent's own.
p <- distributions7:::fold_parts(gaussian1_distrib(), y, theta)
sp <- distrib_gradient(gaussian1_distrib(), y, theta)$mu
sm <- distrib_gradient(gaussian1_distrib(), -y, theta)$mu
all.equal(g$mu, p$w * sp + (1 - p$w) * sm)
#> [1] TRUE
```
