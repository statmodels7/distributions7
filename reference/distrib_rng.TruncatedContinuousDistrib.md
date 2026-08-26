# Truncated Random Number Generator (Continuous)

Draws by inverse transform on the parent, \\Y = Q\\F(L^-) + V Z\\\\ with
\\V\\ uniform on \\(0, 1)\\. The draw is exact and terminates in one
pass however small \\Z\\ is, where rejection sampling from the parent
would need \\1/Z\\ draws on average.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- n:

  Number of observations to generate, a single positive integer.

- theta:

  A named list of the parent's parameters. A component varying by
  observation must have length `n`.

## Value

A numeric vector of `n` random draws, every one inside \\\[L, U\]\\.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_rng()`](https://statmodels7.github.io/distributions7/reference/trunc_rng.md)
for the shared body,
[`distrib_quantile.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedContinuousDistrib.md),
which does the work, and
[`distrib_rng.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedDiscreteDistrib.md)
for the discrete branch.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

set.seed(1)
y <- distrib_rng(tn, 5000, theta)
c(inside = all(y >= -1 & y <= 2), sampled = mean(y),
  truncated = mean(tn, theta), parent = 0.3)
#>    inside   sampled truncated    parent 
#> 1.0000000 0.4099989 0.4159512 0.3000000 

# An interval carrying little mass costs no more than any other.
far <- truncated(gaussian1_distrib(), lower = 4, upper = 5)
set.seed(2)
range(distrib_rng(far, 1000, theta))
#> [1] 4.000226 4.995428
```
