# The Analytic Third and Fourth Derivatives of the Expected Information

Reads a family's kernel at order 3 or 4 on the parameter scale and
carries the result onto the link scale where it is requested, from the
expected information and its derivatives of every lower order.

## Usage

``` r
d34expected_analytic(distrib, y, theta, scale, order, threads, kern)
```

## Arguments

- distrib:

  A distribution object.

- y, theta:

  As the generic takes them.

- scale:

  `"parameter"` or `"link"`.

- order:

  `3L` or `4L`.

- threads:

  The thread count passed to the kernels.

- kern:

  A function of the order returning the kernel's list.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

## See also

[`dexpected_link_higher()`](https://statmodels7.github.io/distributions7/reference/dexpected_link_higher.md)
