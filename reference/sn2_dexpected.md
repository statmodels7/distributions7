# The Centered Skew Normal's Expected Information Derivatives

The parent's expected information and its derivatives, from
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
carried through the map to the centered parametrization by
[`dexpected_chain()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain.md).

## Usage

``` r
sn2_dexpected(distrib, y, theta, order, threads = 1L)
```

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y, theta:

  As the generics take them; the skewness must not be zero.

- order:

  `1L` or `2L`.

- threads:

  The thread count passed to the parent's kernels.

## Value

A named list on the parameter scale, keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).
