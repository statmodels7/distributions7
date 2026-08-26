# Zero-Inflated Random Number Generator

Draws `n` values from the parent, then replaces a Bernoulli(\\\zeta\\)
fraction of them with structural zeros. It consumes the parent's draws
followed by `n` uniforms, so the stream is reproducible under
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) but is not the
parent's stream with a filter: the uniforms come after.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- n:

  The number of draws, a single non-negative whole number.

- theta:

  A named list with the parent's parameters followed by `zi`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `n`.

## Notation

\\f\\ is the parent's mass function, \\\zeta\\ the probability of a
structural zero, \\L_0 = \zeta + (1-\zeta)f(0)\\ the inflated mass at
zero, \\w = (1-\zeta)f(0)/L_0\\ the posterior probability that an
observed zero came from the parent, \\s\\ the parent's score and
\\\ell\\ the log-mass of one observation.

## See also

[`distrib_pdf.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.ZeroInflatedDistrib.md)
for the mass function these are drawn from, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)

set.seed(1)
distrib_rng(d, 10, theta)
#>  [1] 0 0 3 5 2 5 6 4 3 1

# A large sample reproduces the inflated mass at zero, which is well above
# the parent's.
set.seed(1)
big <- distrib_rng(d, 20000, theta)
c(sampled = mean(big == 0), exact = distrib_pdf(d, 0, theta),
  parent = dpois(0, 3))
#>    sampled      exact     parent 
#> 0.28890000 0.28734030 0.04978707 
```
