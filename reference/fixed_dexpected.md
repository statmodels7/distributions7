# The Derivatives of a Fixed Family's Expected Information

Returns the parent family's first (`k = 1`) or second (`k = 2`)
derivatives of the expected information at the full parameter vector,
subset to the free parameters, on the parameter scale. Registered
through
[`register_dexpected()`](https://statmodels7.github.io/distributions7/reference/register_dexpected.md),
which carries them to the link scale with the wrapper's links.

## Usage

``` r
fixed_dexpected(d, y, th, k, t)
```

## Arguments

- d:

  A fixed family, from
  [`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md).

- y:

  The response, read for its length.

- th:

  The free parameters.

- k:

  The order, `1L` or `2L`.

- t:

  The thread count passed to the parent's methods.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
(`k = 1`) or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
(`k = 2`) over the free parameters.
