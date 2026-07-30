# Probability the Parent Puts on a Single Point

\\P(Y = x)\\ under the parent: the pmf for a lattice distribution, the
atom's probability for a mixed one, and zero for an ordinary continuous
distribution.

## Usage

``` r
parent_mass_at(distrib, x, theta)
```

## Arguments

- distrib:

  A truncated distribution object.

- x:

  The point to evaluate at.

- theta:

  A named list of parameters.

## Value

A numeric vector of probabilities.

## Details

This is the one quantity separating the two truncation classes, and
getting it wrong for a *mixed* parent is subtle. It is tempting to
branch on whether the parent is a `discrete_distrib`; that looks right
and is wrong. The cdf of `zero_adjusted(gamma_distrib())` already
includes the point mass at zero, so \\F(0) \neq F(0^-)\\ even though the
object is a `continuous_distrib`. Truncating it from above, with the
atom retained, then drops exactly that mass out of the normalising
constant – and the resulting density integrates to something other than
one while every formula still reads correctly. Asking
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
instead of asking the class cannot make that mistake.

## See also

[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md),
[`trunc_constants`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
