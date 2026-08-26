# Atoms of a Truncated Continuous Distribution

Reports the parent's point masses that survive truncation, rescaled by
\\1/Z\\. An ordinary continuous parent has none and the result is empty;
this matters when the parent is itself MIXED, as
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
of a continuous distribution is.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- theta:

  A named list of the parent's parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `y`, the surviving atom locations, and `p`, their
probabilities under the truncated law. Both are length zero for a parent
with no atoms.

## Details

A mass outside \\\[L, U\]\\ is dropped, along with the rest of the
density there, and one inside is rescaled by the same factor the density
is. Two consumers read the result:
[`expectation.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/expectation.TruncatedContinuousDistrib.md)
adds the masses to its integral, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
adjusts the checks that a point mass invalidates.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for the generic,
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
and
[`distrib_atoms.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.ZeroAdjustedContinuousDistrib.md)
for the parent that declares one.

## Examples

``` r
# An ordinary continuous parent has no atoms.
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)
lengths(distrib_atoms(tn, theta))
#> y p 
#> 0 0 

# A zero-adjusted parent does, and truncation rescales it by 1/Z.
tz <- truncated(zero_adjusted(gamma2_distrib()), upper = 5)
th <- list(mu = 2, sigma2 = 1, za = 0.3)
distrib_atoms(tz, th)
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3021864
#> 
0.3 / distributions7:::trunc_constants(tz, th)$Z
#> [1] 0.3021864
```
