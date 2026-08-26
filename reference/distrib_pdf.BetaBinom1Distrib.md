# Beta-Binomial Probability Mass Function

Computes the beta-binomial mass \$\$P(Y = y) = \binom{n}{y} \dfrac{B(y +
\alpha,\\ n - y + \beta)}{B(\alpha, \beta)}, \qquad \alpha =
\dfrac{\mu}{\sigma}, \quad \beta = \dfrac{1-\mu}{\sigma},\$\$ with \\B\\
the beta function and \\n\\ the object's `size`. Off the support \\\\0,
\dots, n\\\\ the mass is 0, and a non-integer `y` is off the support
too.

The two beta functions are of magnitude
\\(\alpha+\beta)\log(\alpha+\beta)\\ while their difference is of order
one, so writing the mass as that difference loses every digit at a small
\\\sigma\\. The compiled kernel switches to a sum of logarithms instead,
the shifts being integers, and stays exact to the binomial limit.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

- y:

  A numeric vector of counts. A value outside \\\\0, \dots, n\\\\ or not
  an integer gives a mass of 0, or `-Inf` with `log = TRUE`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `mu` must lie in \\(0, 1)\\ and `sigma` be strictly
  positive.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned. Defaults to
  `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(y), length(mu), length(sigma))`, one value per observation.

## See also

[`distrib_cdf.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BetaBinom1Distrib.md)
for the cumulative sum,
[`distrib_gradient.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom1Distrib.md)
for the derivatives of the log-mass,
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic and
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
for the family.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)

# The mass over the whole support sums to one.
m <- distrib_pdf(d, 0:10, th)
round(m, 4)
#>  [1] 0.2645 0.1526 0.1169 0.0965 0.0821 0.0708 0.0612 0.0525 0.0440 0.0350
#> [11] 0.0240
sum(m)
#> [1] 1

# It is the beta-binomial mass at the two implied shapes.
a <- 0.3 / 0.5; b <- 0.7 / 0.5
all.equal(m, choose(10, 0:10) * beta(0:10 + a, 10 - 0:10 + b) / beta(a, b))
#> [1] TRUE

# Off the support, and at a non-integer count, the mass is zero.
distrib_pdf(d, c(-1, 2.5, 11), th)
#> [1] 0 0 0

# As the dispersion goes to zero the family becomes the binomial, and the
# mass stays accurate there where the two beta functions would cancel.
max(abs(distrib_pdf(d, 0:10, list(mu = 0.3, sigma = 1e-6)) -
        dbinom(0:10, 10, 0.3)))
#> [1] 1.334134e-06
```
