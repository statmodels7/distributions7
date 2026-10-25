# The Third and Fourth Derivatives of a Reparametrized Family's Expected Information

Returns the third (`k = 3`) or fourth (`k = 4`) derivatives of the
expected information of a reparametrized family on the parameter scale,
from the parent's closed derivatives and the map's tables. The function
is registered through
[`register_d34expected()`](https://statmodels7.github.io/distributions7/reference/register_d34expected.md),
which carries the derivatives to the link scale with the wrapper's
links.

## Usage

``` r
reparam_d34expected(d, y, th, k, t)
```

## Arguments

- d:

  A reparametrized family, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- y:

  The response, read for its length.

- th:

  The new parameters.

- k:

  The order, `3L` or `4L`.

- t:

  The thread count passed to the parent's methods.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
(`k = 3`) or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)
(`k = 4`).

## Details

The third derivative is the chain of
[`dexpected_chain3()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain3.md)
and the fourth that of
[`dexpected_chain4()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain4.md).
They are exact wherever the parent's derivatives to that order are, and
the map's partials to the fourth order are, the partials of the fifth
order being formed from the fourth by
[`reparam_map_order5()`](https://statmodels7.github.io/distributions7/reference/reparam_map_order5.md).
A family built without `map_derivs` has numerical tables that lose their
leading digits from the third order, and both orders then come from
[`numerical_dexpected_higher()`](https://statmodels7.github.io/distributions7/reference/numerical_dexpected_higher.md),
the single stencil on the second derivative.
