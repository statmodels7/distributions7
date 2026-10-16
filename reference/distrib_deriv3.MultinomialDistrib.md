# Multinomial Third Derivatives

Computes every third derivative of the log-mass in the parameters, in
closed form. Up to a constant in the counts the log-mass is \\\sum_j y_j
\log p_j\\. On the additive log-ratio chart of
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
this equals \\\sum\_{j\<p} y_j\eta_j - N\\\mathrm{LSE}(\eta)\\, with
\\N\\ the total of the row, so every derivative of order two or more is
\\-N\\ times a joint cumulant of the indicators of the categorical law
(see
[`mn_cumulant()`](https://statmodels7.github.io/distributions7/reference/mn_cumulant.md)).
The cumulant is a sum of products of probabilities, computed without
dividing by one, so it stays finite when a probability underflows. On
another chart each term depends on one coordinate of the simplex and the
chain rule reduces to a univariate partition sum per coordinate, with
\$\$f^{(m)}(p) = (-1)^{m-1}(m-1)!\\p^{-m},\$\$ weighted by the counts.

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

  Logical of length 1. When `TRUE` the expectation is returned, in
  closed form (see Details). Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`. The method returns the
  parameter scale and the generic applies the link scale to its result.

- approx:

  Not read: the expectation is in closed form. It is present so that the
  signature matches the generic's.

- nsim:

  Not read. It is present so that the signature matches the generic's.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length \\n\\, keyed and ordered as
`deriv_names(distrib@params, 3)`. At \\p = 3\\ there are four
components, the simplex spending \\p - 1\\ free values. With
`expected = TRUE` every vector is constant.

## Details

With `expected = TRUE` the expectation is exact. The derivatives are
linear in the counts, so the method replaces each count by its
expectation \\n p_j\\. On the additive log-ratio chart this replaces
\\N\\ by \\n\\, and the result equals the observed value at every row on
the support.

## Notation

\\p_j\\ is the probability of category \\j\\, \\y_j\\ its count, \\n\\
the number of trials, \\\eta\\ the simplex's free vector and \\\ell\\
the log-mass of one observation.

## See also

[`distrib_deriv4.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.MultinomialDistrib.md)
for the next order,
[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md)
for the shared engine,
[`mn_cumulant()`](https://statmodels7.github.io/distributions7/reference/mn_cumulant.md)
for the cumulants, and
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

# Against a central difference of the Hessian.
h <- 1e-4
tp <- theta; tp$probs_alr2 <- tp$probs_alr2 + h
tm <- theta; tm$probs_alr2 <- tm$probs_alr2 - h
c(exact = sum(d3[["probs_alr1_probs_alr1_probs_alr2"]]),
  stencil = (sum(distrib_hessian(d, y, tp)[["probs_alr1_probs_alr1"]]) -
             sum(distrib_hessian(d, y, tm)[["probs_alr1_probs_alr1"]])) /
            (2 * h))
#>     exact   stencil 
#> 0.4072999 0.4072999 

# The expectation equals the sum over the finite support.
sup <- mv_support(d, theta)
w <- distrib_pdf(d, sup, theta)
c(expected = distrib_deriv3(d, y, theta, expected = TRUE)[[1]][1],
  sum = sum(w * distrib_deriv3(d, sup, theta)[[1]]))
#>   expected        sum 
#> -0.1914194 -0.1914194 
```
