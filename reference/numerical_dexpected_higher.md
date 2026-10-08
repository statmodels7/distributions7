# The Third or Fourth Derivative of the Expected Information by One Stencil

The default route behind
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md)
and
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md):
one central stencil on the family's analytic second derivative of its
expected information.

## Usage

``` r
numerical_dexpected_higher(distrib, y, theta, scale, order)
```

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  `"parameter"` or `"link"`.

- order:

  `3L` or `4L`.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

## Details

A component of order three, keyed by the pair \\(a,b)\\ and the sorted
triple \\(c,d,e)\\, is the first difference along \\e\\ of the second
derivative in \\(c,d)\\. A component of order four, keyed by the sorted
quadruple \\(c,d,e,f)\\, is the second-order stencil in \\(e,f)\\ of the
second derivative in \\(c,d)\\. Each shifted evaluation is computed once
and shared by every component that reads it: \\2p\\ evaluations at order
three and \\1 + 2p^2\\ at order four.
