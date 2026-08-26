# Multinomial Probability Mass Function

Computes the multinomial mass \$\$P(Y = y) = \dfrac{n!}{\prod\_{j=1}^{p}
y_j!}\prod\_{j=1}^{p} p_j^{y_j},\$\$ one value per row of `y`, on the
weak compositions of \\n\\ into \\p\\ parts. A row with a negative
entry, a non-integer entry, or a sum other than \\n\\ is off the support
and returns 0.

The factorials are formed through
[`base::lgamma()`](https://rdrr.io/r/base/Special.html) and the product
through a matrix multiplication by \\\log p\\, so the whole calculation
runs on the log scale and a mass far below the smallest representable
double still has a finite logarithm.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- y:

  A numeric matrix with one row per observation and \\p\\ columns, each
  row a vector of non-negative integers summing to the object's `size`.
  A single observation may be given as a plain vector and is read as one
  row. A row off the support gives a mass of 0, or `-Inf` with
  `log = TRUE`.

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1. A parameter may not vary by observation here, the
  probability vector being one point of the simplex for the whole
  sample.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned. Defaults to
  `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\ of length `nrow(y)`,
one per observation.

## See also

[`mv_support.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_support.MultinomialDistrib.md)
for the points the mass sits on,
[`distrib_gradient.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MultinomialDistrib.md)
for the derivatives of the log-mass,
[`mv_marginal.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MultinomialDistrib.md)
for a coordinate's binomial marginal, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)

# The mass over the whole support sums to one, exactly.
supp <- mv_support(d, th)
c(points = nrow(supp), mass = sum(distrib_pdf(d, supp, th)))
#> points   mass 
#>     21      1 

# It is stats::dmultinom at the implied probabilities.
pr <- mv_location(d, th) / 5
all.equal(distrib_pdf(d, supp, th),
          apply(supp, 1, function(r) dmultinom(r, prob = pr)))
#> [1] TRUE

# A row that does not sum to the size, or is not integral, is off the
# support.
distrib_pdf(d, rbind(c(1, 1, 1), c(-1, 3, 3), c(1.5, 1.5, 2)), th)
#> [1] 0 0 0
```
