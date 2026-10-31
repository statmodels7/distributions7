# The Fifth Derivative of the Expected Information by One Stencil

Computes the default route behind
[`distrib_d5expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d5expected_hessian.md)
by one tensor-product central stencil of
[`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html)
on the highest derivative of the expected information that the family
computes itself, as
[`d5expected_base()`](https://statmodels7.github.io/distributions7/reference/d5expected_base.md)
gives it.

## Usage

``` r
numerical_d5expected(distrib, y, theta, scale, accuracy = 4L)
```

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  `"parameter"` or `"link"`; on `"link"` the stencil is taken in the
  linear predictors on the link-scale derivatives of the base order.

- accuracy:

  The order of accuracy of the stencils, passed to
  [`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html).

## Value

A named list keyed as
[`d5expected_names()`](https://statmodels7.github.io/distributions7/reference/d5expected_names.md).

## Details

For each pair \\(a,b)\\,
[`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html)
is called with \\\mathbb{E}\[\ell\_{ab}\]\\ as the function and its
derivatives of the base order as the known quantity, and returns every
quintuple. The derivatives of the base order are computed once per node
of the stencil and shared by the pairs. With base 4 and accuracy 4 the
cost is \\4p\\ evaluations of
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md),
for \\p\\ parameters.
