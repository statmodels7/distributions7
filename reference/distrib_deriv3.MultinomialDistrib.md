# Multinomial Third Derivatives

Computes every third derivative of the log-mass in the parameters, in
closed form. Up to a constant in the counts the log-mass is \\\sum_j y_j
\log p_j\\, so each term depends on ONE coordinate of the simplex and
the chain rule collapses to a univariate partition sum per coordinate,
with \$\$f^{(m)}(p) = (-1)^{m-1}(m-1)!\\p^{-m}.\$\$ The combinatorial
factor carries no parameter, so it drops out of every derivative and the
counts enter only as weights.

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
`deriv_names(distrib@params, 3)`. At \\p = 3\\ there are four
components, the simplex spending \\p - 1\\ free values. With
`expected = TRUE` every vector is constant.

## Details

The support is finite, so the expectation of any component is an exact
sum over `mv_support(distrib, theta)` and needs no sampling.
`expected = TRUE` nonetheless routes to
[`mv_expected_higher()`](https://statmodels7.github.io/distributions7/reference/mv_expected_higher.md),
which draws; a caller who wants the exact expectation forms the sum
against
[`distrib_pdf.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MultinomialDistrib.md)
directly.

## Notation

\\p_j\\ is the probability of category \\j\\, \\y_j\\ its count,
\\\eta\\ the simplex's free vector and \\\ell\\ the log-mass of one
observation.

## See also

[`distrib_deriv4.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.MultinomialDistrib.md)
for the next order,
[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md)
for the shared engine,
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
for the exact route to an expectation, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 10)
theta <- list(probs_alr1 = 0.4, probs_alr2 = -0.3)
set.seed(2)
y <- distrib_rng(d, 5, theta)

d3 <- distrib_deriv3(d, y, theta)
vapply(d3, sum, numeric(1))
#> probs_alr1_probs_alr1_probs_alr1 probs_alr1_probs_alr1_probs_alr2 
#>                       -0.9570972                        0.4072999 
#> probs_alr1_probs_alr2_probs_alr2 probs_alr2_probs_alr2_probs_alr2 
#>                        2.8642703                       -4.7842481 

# Against one stencil on the analytic Hessian, component by component.
h <- 1e-4
tp <- theta; tp$probs_alr2 <- tp$probs_alr2 + h
tm <- theta; tm$probs_alr2 <- tm$probs_alr2 - h
c(exact = sum(d3[["probs_alr1_probs_alr1_probs_alr2"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["probs_alr1_probs_alr1"]]) -
             sum(distrib_hessian(d, y, tm)[["probs_alr1_probs_alr1"]])) /
            (2 * h))
#>     exact   stencil 
#> 0.4072999 0.4072999 

# The support is finite, so the mass over it sums to one exactly and an
# expectation can be taken as a sum rather than a sample.
sup <- mv_support(d, theta)
c(states = nrow(sup), mass = sum(distrib_pdf(d, sup, theta)))
#> states   mass 
#>     66      1 
```
