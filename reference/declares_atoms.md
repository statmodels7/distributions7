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
is registered on a class strictly below `distrib` other than the three
passing wrappers, or when one of those wraps a distribution for which
the function returns `TRUE`; `FALSE` when the method comes from the base
class or is absent.

## Details

The wrappers
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
and
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
register
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
on their continuous classes only to report the parent's atoms, so their
registration says nothing about the wrapped family. For these classes
the function is applied to `parent@parent_distrib` instead, and a
wrapper of a family without atoms declares none.

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

# A wrapper declares atoms exactly when the family it wraps does.
distributions7:::declares_atoms(fixed(gaussian1_distrib(), sigma = 1.2))
#> [1] FALSE
distributions7:::declares_atoms(fixed(zero_adjusted(gaussian1_distrib()),
                                      sigma = 1.2))
#> [1] TRUE
```
