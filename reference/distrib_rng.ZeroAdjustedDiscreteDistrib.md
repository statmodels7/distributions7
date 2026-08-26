# Zero-Adjusted Discrete Random Number Generator

Draws zeros with probability \\\pi\\ and otherwise samples from the
ZERO-TRUNCATED parent, so no draw from the positive part is ever zero.
The truncated draw is taken by inverting the parent's distribution
function at a uniform rescaled onto \\(f(0), 1)\\, which needs no
rejection loop and so terminates in bounded time however small the
positive mass is.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- n:

  The number of draws, a single non-negative whole number.

- theta:

  A named list with the parent's parameters followed by `za`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `n`, every value a non-negative support point
of the parent.

## Notation

\\f\\ is the parent's mass function, \\F\\ its distribution function,
\\\pi\\ the probability of a zero, \\s\\ the parent's score and \\\ell\\
the log-mass of one observation.

## See also

[`distrib_pdf.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroAdjustedDiscreteDistrib.md)
for the mass function these are drawn from,
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
for the truncation the positive part uses, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)

set.seed(1)
distrib_rng(d, 10, theta)
#>  [1] 0 0 2 2 0 4 2 4 3 0

# A large sample reproduces the mass at zero, and no positive draw is zero.
set.seed(1)
big <- distrib_rng(d, 20000, theta)
c(sampled_zero = mean(big == 0), parameter = 0.4,
  smallest_positive = min(big[big > 0]))
#>      sampled_zero         parameter smallest_positive 
#>            0.4009            0.4000            1.0000 
```
