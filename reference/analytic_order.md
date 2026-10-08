# The Highest Derivative Order a Family Implements Itself

Returns the highest order, up to `upto`, at which the family registers
its own method: 0 for a family with a density alone, 1 with a score, 2
with a Hessian, 3 with third derivatives. A method inherited from one of
the package's base classes is a numerical fallback and does not count;
one inherited from a family class does.

## Usage

``` r
analytic_order(distrib, upto)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- upto:

  The highest order to look for, 1 to 3.

## Value

An integer from 0 to `upto`.

## See also

[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md),
which differences that order.
