# Truncated Probability Density Function

Evaluates \$\$f_T(y) = \frac{f(y;\theta)}{Z(\theta)} \quad (L \le y \le
U), \qquad f_T(y) = 0 \\ \text{ otherwise},\$\$ with \\Z(\theta) =
F(U;\theta) - F(L;\theta)\\ the retained mass. The division is carried
out on the log scale, so a point far in a tail keeps its precision
instead of underflowing before the division.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of observations. A point outside \\\[L, U\]\\ gives
  `0`, or `-Inf` on the log scale.

- theta:

  A named list of the parent's parameters. Truncation adds none.

- log:

  Logical, default `FALSE`. When `TRUE` the log-density is returned.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of density values, of the length `y` and `theta`
recycle to.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the constructor,
[`trunc_pdf()`](https://statmodels7.github.io/distributions7/reference/trunc_pdf.md)
for the shared body,
[`distrib_pdf.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedDiscreteDistrib.md)
for the discrete branch, and
[`distrib_cdf.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedContinuousDistrib.md).

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

distrib_pdf(tn, c(-2, 0, 1, 3), theta)
#> [1] 0.0000000 0.4118505 0.3584437 0.0000000

# The parent's density lifted by the reciprocal of the retained mass.
Z <- pnorm(2, 0.3, 1.2) - pnorm(-1, 0.3, 1.2)
c(dnorm(0, 0.3, 1.2), dnorm(1, 0.3, 1.2)) / Z
#> [1] 0.4118505 0.3584437

# It integrates to one over the interval; the parent's does not.
integrate(function(y) distrib_pdf(tn, y, theta), -1, 2)$value
#> [1] 1

# Outside the interval the log-density is -Inf, not a small number.
distrib_pdf(tn, c(-2, 3), theta, log = TRUE)
#> [1] -Inf -Inf
```
