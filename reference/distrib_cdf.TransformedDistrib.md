# Transformed Cumulative Distribution Function

Computes \\F_Y(q) = F_X(g^{-1}(q))\\ by evaluating the parent's own
distribution function at the preimage. For a DECREASING transformation
the tails swap: `lower.tail` is inverted before the parent is called,
because \\Y \le q\\ is \\X \ge g^{-1}(q)\\ there.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- q:

  A numeric vector of quantiles on the transformed scale.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, probabilities are \\P(Y
  \le q)\\. For a decreasing transformation this is inverted before
  reaching the parent.

- log.p:

  Logical of length 1. When `TRUE` the logarithm is returned, computed
  by the parent. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, in \\\[0, 1\]\\.

## Details

`lower.tail` and `log.p` are passed THROUGH to the parent and not
applied afterwards, so a parent with an accurate upper tail or an
accurate log-probability keeps that accuracy here.

## Notation

\\g\\ is the transformation and \\F_X\\, \\F_Y\\ the parent's and the
transformed distribution function.

## See also

[`distrib_pdf.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TransformedDistrib.md)
for the density,
[`distrib_quantile.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TransformedDistrib.md),
which inverts this, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)
q <- c(0.5, 1, 3)

distrib_cdf(d, q, theta)
#> [1] 0.0679238 0.2659855 0.7728499
all.equal(distrib_cdf(d, q, theta), plnorm(q, 0.5, 0.8))
#> [1] TRUE

# A decreasing transformation swaps the tails: the reciprocal of a gamma.
ig <- transformation(gamma1_distrib(), inverse_transform())
th2 <- list(mu = 2, phi = 0.3)
ig@transformer@decreasing
#> [1] TRUE
c(transformed = distrib_cdf(ig, 0.5, th2),
  parent_upper = distrib_cdf(gamma1_distrib(), 2, th2, lower.tail = FALSE))
#>  transformed parent_upper 
#>    0.4271256    0.4271256 
```
