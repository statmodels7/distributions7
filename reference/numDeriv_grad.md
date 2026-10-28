# A Central-Difference Gradient Without a Dependency

Computes one central difference per coordinate, which is enough for the
multivariate checks and keeps numDeriv out of the hard dependencies.

## Usage

``` r
numDeriv_grad(f, x, rel_step = NULL)
```

## Arguments

- f:

  A function of a numeric vector, returning one number.

- x:

  The point.

- rel_step:

  Deprecated and unused; the step is
  [`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)'s.

## Value

A numeric vector the length of `x`.

## Details

The nodes, the weights and the step are numericals7's.
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
is not called directly because its `f` maps a vector of points to the
values at those points, while this one maps a whole vector to a single
number.

## See also

[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
