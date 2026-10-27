# Folded Third Derivatives

Computes every third derivative of the folded log-density from the SAME
partition sums that the Hessian uses, one order up. The block ratios are
\$\$\frac{d^B L}{L} = w \frac{d^B f(x)}{f(x)} + (1-w) \frac{d^B
f(-x)}{f(-x)},\$\$ in which each term is a complete Bell polynomial in
the parent's own log-derivatives, and
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
turns them into derivatives of \\\log L\\ by the moment-to-cumulant
relation. A parent with closed forms to third order gives the folded
family closed forms to third order.

## Arguments

- distrib:

  A `FoldedDistrib` object, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- y:

  A numeric vector of observations, non-negative.

- theta:

  A named list of the parent's parameters.

- expected:

  Logical of length 1. When `TRUE` the expectation is approximated by
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`. The generic applies
  the link scale to the value that the method returns.

- method:

  `NULL` (the default) or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  It is read only when `expected` is `TRUE`: `NULL` selects
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md),
  and an object selects its own route. The generic applies an object
  before dispatching on the family, so that the object is used also by a
  family that has a closed form. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, keyed and ordered as
`deriv_names(distrib@params, 3)`.

## Details

With `expected = TRUE` the expectation goes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
which reads `method`; there is no closed-form expectation for a fold in
general.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one, \\B\\ a multiset of
parameter names and \\w\\ the weight of the positive preimage.

## See also

[`distrib_deriv4.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.FoldedDistrib.md)
for the next order,
[`distrib_hessian.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.FoldedDistrib.md)
for the second,
[`fold_deriv_k()`](https://statmodels7.github.io/distributions7/reference/fold_deriv_k.md)
for the shared body, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)
set.seed(2)
y <- distrib_rng(d, 30, theta)

d3 <- distrib_deriv3(d, y, theta)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# The result is compared with one stencil on the analytic Hessian.
h <- 1e-4
tp <- theta; tp$sigma <- tp$sigma + h
tm <- theta; tm$sigma <- tm$sigma - h
c(exact = sum(d3[["mu_mu_sigma"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["mu_mu"]]) -
             sum(distrib_hessian(d, y, tm)[["mu_mu"]])) / (2 * h))
#>     exact   stencil 
#> -1.241218 -1.241217 
```
