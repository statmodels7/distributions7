# Does a Distribution Declare Atoms

Answers whether a distribution registers
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for itself, which is how
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
decides that a parent carries a point mass and must be rejected. The
question is asked of the CLASS the method is registered on, not of a
value: a distribution with an atom declares one at every parameter
setting, and a value taken at one setting could be empty by accident.

## Usage

``` r
declares_atoms(parent)
```

## Arguments

- parent:

  A `distrib` object.

## Value

`TRUE` when
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
is registered on a class strictly below `distrib`, `FALSE` when the
method comes from the base class or is absent.

## Details

The argument is named `parent` deliberately. The base class of this
package is called `distrib`, and an argument of that name would shadow
it: the comparison meant for the base class would then be against the
object. That defect has been met before in this package and is avoided
by naming.

## See also

[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for the generic,
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md),
which consults this, and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
which produces a parent it rejects.

## Examples

``` r
# A plain gaussian declares none.
distributions7:::declares_atoms(gaussian1_distrib())
#> [1] FALSE

# A zero-adjusted continuous parent is a mixed distribution and does.
distributions7:::declares_atoms(zero_adjusted(gaussian1_distrib()))
#> [1] TRUE

# Which is why folded() rejects the second by name.
try(folded(zero_adjusted(gaussian1_distrib())))
#> Error : 'zero-adjusted gaussian1' carries an atom, and folding would misplace it: zero is its own
#>   preimage while every other point has two, so an atom at zero would be
#>   counted twice and one elsewhere moved onto its reflection.
```
