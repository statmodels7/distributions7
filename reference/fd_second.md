# A Second Derivative From One Stencil

The mixed or repeated second derivative of a scalar function, from a
single stencil rather than from nested first differences.

## Usage

``` r
fd_second(f, x, k, l, h_rel = NULL)
```

## Arguments

- f:

  A function of a numeric vector, returning one number.

- x:

  The point.

- k, l:

  The coordinates to differentiate in.

- h_rel:

  Deprecated and unused; the step is
  [`fd_step`](https://statmodels7.github.io/numericals7/reference/fd_step.html)'s
  at order two.

## Value

A single number.

## Details

The nodes, the weights and the step are numericals7's. Where the two
coordinates differ the stencil is the product of two first-order
factors, which is one stencil in two variables and not a difference of a
difference: nesting is forbidden along ONE variable, and is what a mixed
derivative is along two.

## See also

[`fd_weights`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
