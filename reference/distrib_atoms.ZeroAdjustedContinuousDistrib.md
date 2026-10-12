# Atoms of a Zero-Adjusted Continuous Distribution

Declares the single point mass at zero, with probability \\\pi\\. This
declaration marks the object as a MIXED distribution: its density
integrates to \\1 - \pi\\, and the remaining probability \\\pi\\ is the
mass of the atom. A function that ignores the atom finds a total
probability of \\1 - \pi\\.

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

The declaration is read by
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
and
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
changes the checks that fail on a mixed law: the density integrates to
\\1-\pi\\ plus the atom, the quantile round trip cannot close inside the
jump, and a central difference in \\y\\ across the atom has no
derivative to converge to.
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
keeps the atoms of the parent that lie inside the truncation interval.
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
rejects a parent that declares an atom, because zero is its own preimage
under the absolute value and every other point has two.
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
does not read the declaration; its method for this class splits the
integral at zero.

The probability is read from `theta`, so the returned value changes with
the value of `za`.

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

# The declaration is what folded() rejects.
try(folded(d))
#> Error : 'zero-adjusted gaussian1' carries an atom, and folding would misplace it: zero is its own
#>   preimage while every other point has two, so an atom at zero would be
#>   counted twice and one elsewhere moved onto its reflection.
```
