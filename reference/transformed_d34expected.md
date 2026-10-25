# The Third and Fourth Derivatives of a Transformed Family's Expected Information

Returns the parent family's third (`k = 3`) or fourth (`k = 4`)
derivatives of the expected information on the parameter scale. The
function is registered through
[`register_d34expected()`](https://statmodels7.github.io/distributions7/reference/register_d34expected.md),
which carries the derivatives to the link scale with the wrapper's
links.

## Usage

``` r
transformed_d34expected(d, y, th, k, t)
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

A transformation of the response leaves the expected information
unchanged as a function of the parameters: the log-density of \\Y =
g(X)\\ differs from the parent's at \\g^{-1}(y)\\ by \\\log\|d
g^{-1}/dy\|\\, which does not depend on the parameters, and \\Y\\ and
\\X\\ have the same distribution of the score. The expected Hessian
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
of the wrapper is accordingly the parent's, and so are its derivatives
in the parameters at every order. The response is carried to the
parent's scale by the inverse transformation, as
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
does, and is read only for its length.
