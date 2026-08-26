# Beta-Binomial Mass Function in Its Shapes

Computes the beta-binomial mass in the two shapes, \$\$P(Y = y) =
\binom{n}{y}\dfrac{B(y+\alpha,\\ n-y+\beta)}{B(\alpha,\beta)},\$\$ with
\\B\\ the beta function and \\n\\ the object's `size`. The support is
tested first, so a count outside \\\\0, \dots, n\\\\ or not an integer
returns a mass of 0 rather than reaching
[`base::lchoose()`](https://rdrr.io/r/base/Special.html), which would
warn on a non-integer argument.

The evaluation route is chosen by the concentration; see
[`betabinom_log_mass()`](https://statmodels7.github.io/distributions7/reference/betabinom_log_mass.md),
which this calls.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

- y:

  A numeric vector of counts. A value outside \\\\0, \dots, n\\\\ or not
  an integer gives a mass of 0, or `-Inf` with `log = TRUE`.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned. Defaults to
  `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`length(y)`, one value per observation.

## See also

[`betabinom_log_mass()`](https://statmodels7.github.io/distributions7/reference/betabinom_log_mass.md)
for the two evaluation routes,
[`distrib_gradient.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom2Distrib.md)
for the derivatives of the log-mass,
[`distrib_pdf.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaBinom1Distrib.md)
for the same quantity in the mean and dispersion, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
th <- list(alpha = 2, beta = 3)

# The mass over the whole support sums to one.
m <- distrib_pdf(d, 0:10, th)
round(m, 4)
#>  [1] 0.0659 0.1099 0.1349 0.1439 0.1399 0.1259 0.1049 0.0799 0.0539 0.0300
#> [11] 0.0110
sum(m)
#> [1] 1

# It is the beta-binomial mass written out.
all.equal(m, choose(10, 0:10) * beta(0:10 + 2, 10 - 0:10 + 3) / beta(2, 3))
#> [1] TRUE

# Off the support, and at a non-integer count, the mass is zero.
distrib_pdf(d, c(-1, 2.5, 11), th)
#> [1] 0 0 0

# The same law as betabinom1 at mu = alpha / (alpha + beta) and
# sigma = 1 / (alpha + beta).
all.equal(m, distrib_pdf(betabinom1_distrib(size = 10), 0:10,
                         list(mu = 0.4, sigma = 0.2)))
#> [1] TRUE
```
