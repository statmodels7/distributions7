# Density of a Truncated Distribution

Evaluates the parent's density divided by the retained mass \\Z\\, and
zero outside \\\[L, U\]\\. This is one of the shared method bodies:
truncation treats the two kinds of parent identically once
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
has resolved the one place they differ, so the body is written once and
registered on both classes.

## Usage

``` r
trunc_pdf(distrib, y, theta, log = FALSE, ...)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- y:

  A numeric vector of observations. A point outside \\\[L, U\]\\ gives
  `0`, or `-Inf` on the log scale.

- theta:

  A named list of the parent's parameters.

- log:

  Logical, default `FALSE`. When `TRUE` the log-density is returned.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of the length `y` and `theta` recycle to.

## Details

The division is done on the LOG scale and the outside is set to
\\-\infty\\ there, so a point far in a tail keeps its precision instead
of underflowing to zero before the division.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`distrib_pdf.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedContinuousDistrib.md)
and
[`distrib_pdf.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedDiscreteDistrib.md),
the two registrations, and
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
for \\Z\\.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

distrib_pdf(tn, c(-2, 0, 1, 3), theta)
#> [1] 0.0000000 0.4118505 0.3584437 0.0000000

# The parent's density divided by the retained mass.
Z <- distributions7:::trunc_constants(tn, theta)$Z
dnorm(c(0, 1), 0.3, 1.2) / Z
#> [1] 0.4118505 0.3584437

# Outside the interval, zero, and -Inf on the log scale.
distrib_pdf(tn, c(-2, 3), theta, log = TRUE)
#> [1] -Inf -Inf

# It integrates to one over the interval, which the parent's does not.
integrate(function(y) distrib_pdf(tn, y, theta), -1, 2)$value
#> [1] 1
```
