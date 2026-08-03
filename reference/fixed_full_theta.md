# Splice the Fixed Values Back Into a Full Parameter List

Combines the wrapper's free `theta` with its fixed values into the full
parameter list of the parent, in the parent's order.

## Usage

``` r
fixed_full_theta(distrib, theta)
```

## Arguments

- distrib:

  A fixed-parameter wrapper object.

- theta:

  A named list or vector of the free parameters.

## Value

A named list covering every parameter of the parent.

## Details

`theta` is aligned against the wrapper first, so the function is safe to
call both from generic-dispatched methods, whose `theta` is already
aligned, and from delegating methods such as
[`mean()`](https://rdrr.io/r/base/mean.html), whose `theta` arrives as
the caller wrote it. Free values may be vectors – the wrapper is as
vectorised in `theta` as its parent – while the fixed values are scalars
by construction.

## See also

[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.md)
