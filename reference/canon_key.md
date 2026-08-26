# Canonical Component Name of a Block

Names a block of parameters the way
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
does: in the order the distribution declares them, joined by an
underscore.

## Usage

``` r
canon_key(block, params)
```

## Arguments

- block:

  A character vector of parameter names.

- params:

  The distribution's parameter names, in declaration order.

## Value

A single string.
