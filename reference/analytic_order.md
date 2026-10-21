# The Highest Derivative Order That a Family Implements Itself

Returns the highest order, up to `upto`, at which the family registers
its own method: 0 for a family that has only a density, 1 with a score,
2 with a Hessian, 3 with third derivatives, 4 with fourth derivatives. A
method inherited from one of the package's base classes is a numerical
fallback and does not count; one inherited from a family class does.

## Usage

``` r
analytic_order(distrib, upto)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- upto:

  The highest order to be examined, 1 to 4.

## Value

An integer from 0 to `upto`.

## See also

[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md),
which differences that order.
