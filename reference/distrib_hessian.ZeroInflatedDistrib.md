# Zero-Inflated Observed Hessian

Computes the second derivatives of the zero-inflated log-mass in closed
form. At a POSITIVE observation the parent block is the parent's own
Hessian, the mixed block is zero and the \\\zeta\\ block is
\\-1/(1-\zeta)^2\\. At a ZERO every block picks up the mixture's
correction: \$\$\ell^{(\theta_i\theta_j)} = w\\H\_{ij}(0) +
w(1-w)\\s_i(0)s_j(0), \qquad \ell^{(\theta_i\zeta)} =
-\frac{f(0)\\s_i(0)}{L_0^2}, \qquad \ell^{(\zeta\zeta)} = -\frac{(1 -
f(0))^2}{L_0^2},\$\$ with \\w = (1-\zeta)f(0)/L_0\\ the posterior weight
of the parent component.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with the parent's parameters followed by `zi`.

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

The parent block is the ordinary two-component mixture curvature: the
weighted Hessian plus the variance of the score across the two
components, \\w(1-w)\\ being that variance's weight when one component's
score is zero.

The mixed block collapses. Written out it is \\-f'(0)/L_0 -
(1-\zeta)f'(0)(1 - f(0))/L_0^2\\, and the bracket is exactly one, so the
whole thing is \\-f'(0)/L_0^2\\.

## Notation

\\f\\ is the parent's mass function, \\\zeta\\ the probability of a
structural zero, \\L_0 = \zeta + (1-\zeta)f(0)\\ the inflated mass at
zero, \\w = (1-\zeta)f(0)/L_0\\ the posterior probability that an
observed zero came from the parent, \\s\\ the parent's score, \\H\\ its
observed Hessian and \\\ell\\ the log-mass of one observation.

## See also

[`distrib_gradient.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroInflatedDistrib.md)
for the first order,
[`distrib_expected_hessian.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroInflatedDistrib.md)
for the expectation, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)
set.seed(2)
y <- distrib_rng(d, 200, theta)

H <- distrib_hessian(d, y, theta)
vapply(H, sum, numeric(1))
#>      mu_mu      zi_zi      mu_zi 
#>  -37.75413 -987.45601   41.60759 

# Against a numerical Hessian of the log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
Hn <- numDeriv::hessian(ll, unlist(theta))
ref <- vapply(distributions7:::hess_pairs(d@params),
              function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
              numeric(1))
max(abs(vapply(H, sum, numeric(1)) - ref))
#> [1] 4.271214e-10

# At a positive observation the mixed block is exactly zero and the zi
# block is a constant; at a zero neither is.
rbind(positive = c(mu_zi = H$mu_zi[y > 0][1], zi_zi = H$zi_zi[y > 0][1]),
      at_zero = c(mu_zi = H$mu_zi[y == 0][1], zi_zi = H$zi_zi[y == 0][1]))
#>              mu_zi      zi_zi
#> positive 0.0000000  -1.777778
#> at_zero  0.6030085 -10.935755
-1 / 0.75^2
#> [1] -1.777778
```
