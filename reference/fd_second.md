# A Second Derivative From One Stencil

The mixed or repeated second derivative of a scalar function, from a
single stencil rather than from nested first differences.

## Usage

``` r
fd_second(f, x, k, l, h_rel = .Machine$double.eps^(1/4))
```

## Arguments

- f:

  A function of a numeric vector, returning one number.

- x:

  The point.

- k, l:

  The coordinates to differentiate in.

- h_rel:

  The relative step.

## Value

A single number.
