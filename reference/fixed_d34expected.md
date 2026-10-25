# The Third and Fourth Derivatives of a Fixed Family's Expected Information

Returns the parent family's third (`k = 3`) or fourth (`k = 4`)
derivatives of the expected information at the full parameter vector,
subset to the free parameters, on the parameter scale. The function is
registered through
[`register_d34expected()`](https://statmodels7.github.io/distributions7/reference/register_d34expected.md),
which carries the derivatives to the link scale with the wrapper's
links.

## Usage

``` r
fixed_d34expected(d, y, th, k, t)
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

  The order, `3L` or `4L`.

- t:

  The thread count passed to the parent's methods.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
(`k = 3`) or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)
(`k = 4`) over the free parameters.

## Details

The expected information of the wrapper is the parent's submatrix over
the free parameters, evaluated at the spliced vector, and the fixed
values do not move. A derivative in free parameters of an entry in free
parameters is therefore the parent's component with every index free.
The free set keeps the parent's order, so the subset by
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)
of the free parameters pairs each component with the same indices as in
the parent. The parent's own route is used, closed form where the parent
has one and the stencil of
[`numerical_dexpected_higher()`](https://statmodels7.github.io/distributions7/reference/numerical_dexpected_higher.md)
otherwise, so the accuracy of each component is the parent's. A fully
fixed family has no free parameter and returns an empty list.
