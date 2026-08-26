# Folded Fourth Derivatives

Computes every fourth derivative of the folded log-density from the same
partition sums as
[`distrib_deriv3.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.FoldedDistrib.md),
one order up: the block ratios \\d^B L/L = w\\(d^B f(x)/f(x)) +
(1-w)\\(d^B f(-x)/f(-x))\\ handed to
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md).
A parent with closed forms to fourth order gives the folded family
closed forms to fourth order.

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

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- approx:

  The strategy
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  uses when `expected = TRUE`. Ignored otherwise.

- nsim:

  The Monte Carlo sample size when `approx = "mc"`. Defaults to `10000`.
  Ignored otherwise.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, keyed and ordered as
`deriv_names(distrib@params, 4)`.

## Details

With `expected = TRUE` the expectation goes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
which reads `approx` and `nsim`.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one, \\B\\ a multiset of
parameter names and \\w\\ the weight of the positive preimage.

## See also

[`distrib_deriv3.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.FoldedDistrib.md)
for the order below,
[`fold_deriv_k()`](https://statmodels7.github.io/distributions7/reference/fold_deriv_k.md)
for the shared body, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)
set.seed(2)
y <- distrib_rng(d, 30, theta)

d4 <- distrib_deriv4(d, y, theta)
length(d4)
#> [1] 5

# Against one stencil on the analytic third order.
h <- 1e-4
tp <- theta; tp$sigma <- tp$sigma + h
tm <- theta; tm$sigma <- tm$sigma - h
c(exact = sum(d4[["mu_mu_sigma_sigma"]]),
  stencil = (sum(distrib_deriv3(d, y, tp)[["mu_mu_sigma"]]) -
             sum(distrib_deriv3(d, y, tm)[["mu_mu_sigma"]])) / (2 * h))
#>     exact   stencil 
#> -68.57366 -68.57366 
```
