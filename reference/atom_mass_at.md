# Mass of an Atom at One Point, by Observation

Returns the probability that a distribution with atoms places on the
single point `x`, one value per observation.
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
reports the masses at the first observation's parameters, so a parameter
that varies by observation is read one row at a time.

## Usage

``` r
atom_mass_at(parent, x, theta)
```

## Arguments

- parent:

  A distribution object that may declare atoms.

- x:

  A single number, the point.

- theta:

  A named list of the parent's parameters. A component may vary by
  observation.

## Value

A numeric vector with one mass per observation (a single number for
scalar parameters), zero where `x` is not an atom.

## See also

[`parent_mass_at()`](https://statmodels7.github.io/distributions7/reference/parent_mass_at.md),
which calls it for a continuous parent, and
[`expectation.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/expectation.TruncatedContinuousDistrib.md).

## Examples

``` r
za <- zero_adjusted(gamma2_distrib())
distributions7:::atom_mass_at(za, 0, list(mu = c(1, 4), sigma2 = 1,
                                          za = c(0.3, 0.6)))
#> [1] 0.3 0.6
```
