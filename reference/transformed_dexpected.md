# The Derivatives of a Transformed Family's Expected Information

Returns the parent family's first (`k = 1`) or second (`k = 2`)
derivatives of the expected information on the parameter scale. A
transformation of the response leaves the expected information unchanged
as a function of the parameters, so its derivatives are the parent's.
Registered through
[`register_dexpected()`](https://statmodels7.github.io/distributions7/reference/register_dexpected.md).

## Usage

``` r
transformed_dexpected(d, y, th, k, t)
```

## Arguments

- d:

  A transformed family, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- y:

  The response, read for its length.

- th:

  The parameters.

- k:

  The order, `1L` or `2L`.

- t:

  The thread count passed to the parent's methods.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
(`k = 1`) or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
(`k = 2`).
