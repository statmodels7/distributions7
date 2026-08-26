# Multinomial Fourth Derivatives

Computes every fourth derivative of the log-mass in the parameters, in
closed form, by the same construction as
[`distrib_deriv3.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MultinomialDistrib.md)
one order up: a univariate partition sum per coordinate of the simplex,
with \\f^{(m)}(p) = (-1)^{m-1}(m-1)!\\p^{-m}\\ and the counts as
weights.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- y:

  A numeric matrix with one row per observation and one column per
  category, each row a vector of counts summing to `distrib@size`.

- theta:

  A named list of parameters, each component a single number:
  `probs_alr1`, ..., `probs_alr(p-1)`.

- expected:

  Logical of length 1. When `TRUE` the expectation is returned, by
  SAMPLING through
  [`mv_expected_higher()`](https://statmodels7.github.io/distributions7/reference/mv_expected_higher.md).
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- approx:

  Ignored: sampling is the only multivariate route to an expectation.
  Present so that the signature matches the generic's.

- nsim:

  The number of draws used when `expected = TRUE`. Defaults to `10000`.
  Ignored otherwise.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 4)`. At \\p = 3\\ there are five
components. With `expected = TRUE` every vector is constant.

## Details

The license for this order is that the SAME assembly run at orders one
and two reproduces the hand-written score and information, which are
derived separately and are already under
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md).

## Notation

\\p_j\\ is the probability of category \\j\\, \\y_j\\ its count,
\\\eta\\ the simplex's free vector and \\\ell\\ the log-mass of one
observation.

## See also

[`distrib_deriv3.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MultinomialDistrib.md)
for the order below,
[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md)
for the shared engine, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 10)
theta <- list(probs_alr1 = 0.4, probs_alr2 = -0.3)
set.seed(2)
y <- distrib_rng(d, 5, theta)

d4 <- distrib_deriv4(d, y, theta)
names(d4)
#> [1] "probs_alr1_probs_alr1_probs_alr1_probs_alr1"
#> [2] "probs_alr1_probs_alr1_probs_alr1_probs_alr2"
#> [3] "probs_alr1_probs_alr1_probs_alr2_probs_alr2"
#> [4] "probs_alr1_probs_alr2_probs_alr2_probs_alr2"
#> [5] "probs_alr2_probs_alr2_probs_alr2_probs_alr2"

# Against one stencil on the analytic third order.
h <- 1e-4
tp <- theta; tp$probs_alr2 <- tp$probs_alr2 + h
tm <- theta; tm$probs_alr2 <- tm$probs_alr2 - h
c(exact = sum(d4[["probs_alr1_probs_alr1_probs_alr2_probs_alr2"]]),
  stencil = (sum(distrib_deriv3(d, y, tp)[["probs_alr1_probs_alr1_probs_alr2"]]) -
             sum(distrib_deriv3(d, y, tm)[["probs_alr1_probs_alr1_probs_alr2"]])) /
            (2 * h))
#>    exact  stencil 
#> 1.339099 1.339099 
```
