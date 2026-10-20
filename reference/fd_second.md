# A Second Derivative From One Stencil

Computes the mixed or repeated second derivative of a scalar function
from a single stencil, without nested first differences.

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

  The coordinates with respect to which `f` is differentiated.

- h_rel:

  Deprecated and unused; the step is
  [`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)'s
  at order two.

## Value

A single number.

## Details

The nodes, the weights and the step are those of numericals7. Where the
two coordinates differ the stencil is the product of two first-order
factors, which is one stencil in two variables and not a difference of a
difference. Nested differences are not used along one variable; a mixed
derivative differentiates along two different variables.

## See also

[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
