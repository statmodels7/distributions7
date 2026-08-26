# Generalized Gamma Distribution Function

Computes \\F(q) = P(d/p,\\ (q/a)^{p})\\, the regularized lower
incomplete gamma function. The identity is exact: \\u = (Y/a)^p\\ is
Gamma with shape \\d/p\\ and unit rate, so the distribution function is
that Gamma's evaluated at \\(q/a)^p\\ and the whole computation is one
`pgamma` call.

`lower.tail` and `log.p` are passed straight to `pgamma`, so the
survival function and the log-probability keep the accuracy R's own
routine gives them in the far tail.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- q:

  A numeric vector of quantiles. A negative value is clamped to zero
  before the transformation, so the value there is 0.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `q`. All three must be strictly
  positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, the value is \\P(Y \le
  q)\\; when `FALSE` it is \\P(Y \> q)\\, computed by `pgamma` rather
  than as one minus the other.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
with `log.p = TRUE`, of the length of the recycled inputs.

## Notation

\\a\\ is the scale, \\d\\ and \\p\\ the shapes, and \\P(s, x) =
\gamma(s,x)/\Gamma(s)\\ the regularized lower incomplete gamma function,
which is `pgamma(x, shape = s)` in R.

## See also

[`distrib_quantile.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.GenGamma1Distrib.md),
which inverts this,
[`distrib_pdf.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.GenGamma1Distrib.md)
for the density, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- gengamma1_distrib()
q <- c(0.5, 1.5, 4)
th <- list(a = 2, d = 3, p = 1.5)

# The identity written out.
all.equal(distrib_cdf(d, q, th), pgamma((q / 2)^1.5, shape = 3 / 1.5))
#> [1] TRUE

# Against a direct quadrature of the density.
rbind(gamma_identity = distrib_cdf(d, q, th),
      quadrature = vapply(q, function(u)
        integrate(function(v) distrib_pdf(d, v, th), 0, u)$value, 0))
#>                       [,1]      [,2]     [,3]
#> gamma_identity 0.007190985 0.1384613 0.773718
#> quadrature     0.007190985 0.1384613 0.773718

# The upper tail keeps its digits where one minus the lower one would not.
c(upper = distrib_cdf(d, 20, th, lower.tail = FALSE),
  one_minus_lower = 1 - distrib_cdf(d, 20, th))
#>           upper one_minus_lower 
#>    6.024535e-13    6.024070e-13 
```
