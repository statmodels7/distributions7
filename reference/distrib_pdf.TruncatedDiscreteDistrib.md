# Truncated Probability Mass Function

Evaluates \$\$P_T(Y = y) = \frac{f(y;\theta)}{Z(\theta)} \quad (L \le y
\le U), \qquad 0 \\ \text{ otherwise},\$\$ with \\Z(\theta) =
F(U;\theta) - F(L;\theta) + f(L;\theta)\\. The mass at the lower
endpoint is added back because that endpoint is INCLUDED in the
truncated support, so `truncated(poisson_distrib(), lower = 1)` is the
zero-truncated Poisson and keeps the mass at one.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of observations. A point outside \\\[L, U\]\\ gives
  `0`, or `-Inf` on the log scale.

- theta:

  A named list of the parent's parameters. Truncation adds none.

- log:

  Logical, default `FALSE`. When `TRUE` the log-probability is returned.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, of the length `y` and `theta` recycle
to.

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
[`distrib_pdf.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TruncatedContinuousDistrib.md)
for the continuous branch, and
[`distrib_cdf.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedDiscreteDistrib.md).

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
theta <- list(mu = 2)

# Zero is gone; one is not.
distrib_pdf(ztp, 0:4, theta)
#> [1] 0.0000000 0.3130353 0.3130353 0.2086902 0.1043451
dpois(1:4, 2) / (1 - dpois(0, 2))
#> [1] 0.3130353 0.3130353 0.2086902 0.1043451

# The retained mass sums to one over the truncated support.
sum(distrib_pdf(ztp, 1:400, theta))
#> [1] 1
```
