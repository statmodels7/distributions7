# Does a Distribution Declare Atoms

Answers whether a distribution registers
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for itself rather than inheriting the base class's empty answer.

## Usage

``` r
declares_atoms(parent)
```

## Arguments

- parent:

  A distribution object.

## Value

A single logical.

## Details

The question is asked of the class rather than of a value, since a
constructor has no `theta` to evaluate the generic at, and the answer is
structural in any case. The class a method was registered on is what
settles it: [`identical()`](https://rdrr.io/r/base/identical.html) on
the method object would not, S7 wrapping it.

The argument is named `parent` deliberately. The base class of this
package is called `distrib`, so an argument of that name would shadow it
and the comparison below would test the class against the distribution
object instead – answering `TRUE` for every family, which is exactly
what it did before this helper existed.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md),
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
