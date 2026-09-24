# The Analytic Derivatives of the Expected Information, on Either Scale

Reads a family's compiled kernel at order 1 or 2 on the parameter scale
and carries the result onto the link scale where it is asked for.

## Usage

``` r
dexpected_analytic(distrib, y, theta, scale, order, threads, kern)
```

## Arguments

- distrib:

  A distribution object.

- y, theta:

  As the generic takes them.

- scale:

  `"parameter"` or `"link"`.

- order:

  `1L` for
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
  `2L` for
  [`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md).

- threads:

  The thread count passed to the kernel and to the family's expected
  information.

- kern:

  A function of the order returning the kernel's named list on the
  parameter scale. A kernel asked for order 2 may return the order-1
  components beside the order-2 ones, and where it does they are read
  from that one call rather than from a second: a route that sums over
  the support pays for the family's derivatives once.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## See also

[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md)
