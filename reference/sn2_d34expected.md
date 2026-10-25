# The Third or Fourth Derivative of the Centered Skew Normal's Expected Information

Computes the components by
[`loc_scale_d34expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_d34expected.md)
where \\\|\gamma_1\|\\ is at least the threshold of the series that
orders zero to two use near zero, and by one stencil on the second
derivative of the expected information, as the default method does, at
the observations below it.

## Usage

``` r
sn2_d34expected(d, y, th, k, t)
```

## Arguments

- d:

  A `skewnormal2` distribution object.

- y:

  A numeric vector of observations, read for its length.

- th:

  A named list of parameters.

- k:

  `3L` or `4L`.

- t:

  The thread count.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md),
on the parameter scale.

## Details

Near \\\gamma_1 = 0\\ the derivatives of the log-density in \\\gamma_1\\
diverge in the representation that the quadrature uses, while the second
derivative of the expected information comes from a series that is exact
there.
