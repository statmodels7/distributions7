# Atoms of a Zero-Adjusted Continuous Distribution

Declares the single point mass at zero, with probability \\\pi\\. This
declaration is what marks the object a MIXED distribution: its density
integrates to \\1 - \pi\\, so a consumer that does not know about the
atom will find the missing \\\pi\\ and call it an error.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- theta:

  A named list with the parent's parameters followed by `za`. Only `za`
  is read, and only its first element.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `y`, the numeric vector `0`, and `p`, the probability
at it.

## Details

Three things read it.
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
adjusts four of its thirteen checks, since the density integrates to
\\1-\pi\\, the quantile round trip cannot close inside the jump, and a
central difference in \\y\\ across the atom has no derivative to
converge to.
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
splits its integral at the atom. And
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
REFUSES a parent that declares one, zero being its own preimage under
the absolute value while every other point has two.

The question is asked of the OBJECT, so the returned probability moves
with `theta`.

## Notation

\\\pi\\ is the probability of the atom at zero.

## See also

[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for the generic,
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
and
[`expectation.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/expectation.ZeroAdjustedContinuousDistrib.md)
for the two consumers, and
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md),
which rejects a parent that declares one.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)

distrib_atoms(d, theta)
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3
#> 

# A plain continuous parent declares none.
distrib_atoms(gaussian1_distrib(), theta[c("mu", "sigma")])
#> $y
#> numeric(0)
#> 
#> $p
#> numeric(0)
#> 

# The declaration is what folded() refuses on.
try(folded(d))
#> Error : 'zero-adjusted gaussian1' carries an atom, and folding would misplace it: zero is its own
#>   preimage while every other point has two, so an atom at zero would be
#>   counted twice and one elsewhere moved onto its reflection.
```
