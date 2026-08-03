# A Central-Difference Gradient Without a Dependency

One central difference per coordinate, which is all the multivariate
checks need and which keeps numDeriv a suggestion rather than a
requirement.

## Usage

``` r
numDeriv_grad(f, x, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- f:

  A function of a numeric vector, returning one number.

- x:

  The point.

- h_rel:

  The relative step.

## Value

A numeric vector the length of `x`.
