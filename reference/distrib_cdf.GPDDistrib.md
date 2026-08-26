# Generalized Pareto Distribution Function

Computes the generalized Pareto distribution function \$\$F(q; \sigma,
\xi) = 1 - \left(1 + \dfrac{\xi q}{\sigma}\right)^{-1/\xi},\$\$ with
\\1 - e^{-q/\sigma}\\ at \\\xi = 0\\. The survival function is formed
first and the tail and the logarithm applied to it, so the far tail
keeps its digits.

The value is 0 below \\q = 0\\ and 1 at and beyond
[`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md),
and is clamped to \\\[0, 1\]\\ so that rounding cannot return a
probability outside the unit interval.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `q`.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, the value is \\P(Y \le
  q)\\; when `FALSE` it is the survival function, which is the quantity
  actually computed.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
with `log.p = TRUE`, of the length of the recycled inputs.

## Details

The branch on \\\|\xi\| \< 10^{-8}\\ is taken by indexing, and that is
not a stylistic choice. [`ifelse()`](https://rdrr.io/r/base/ifelse.html)
returns a result the length of its **test**, so with it a scalar shape
beside a vector of quantiles would collapse the answer to one number.
The test is recycled to the answer's length first.

## Notation

\\\sigma \> 0\\ is the scale and \\\xi\\ the shape.

## See also

[`distrib_quantile.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GPDDistrib.md),
which inverts this in closed form,
[`distrib_pdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GPDDistrib.md)
for the density, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
q <- c(0.2, 1, 4)
th <- list(sigma = 1.5, xi = 0.3)

# The formula written out.
all.equal(distrib_cdf(d, q, th), 1 - (1 + 0.3 * q / 1.5)^(-1 / 0.3))
#> [1] TRUE

# A scalar shape with a vector of quantiles returns one value per
# quantile, which an ifelse() branch would not.
length(distrib_cdf(d, seq(0, 5, by = 1), th))
#> [1] 6

# Shape zero is the exponential.
all.equal(distrib_cdf(d, q, list(sigma = 1.5, xi = 0)),
          pexp(q, rate = 1 / 1.5))
#> [1] TRUE

# A negative shape reaches one at the endpoint and stays there.
distrib_cdf(d, c(4, 5, 6), list(sigma = 2, xi = -0.4))
#> [1] 0.9821115 1.0000000 1.0000000
```
