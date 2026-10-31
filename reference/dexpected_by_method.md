# Derivatives of the Expected Information Under an Explicit Method

Returns the derivative of order `m` of the matrix that an
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
defines, with respect to the parameters, on the parameter scale or the
link scale. It is the route of
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md),
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md),
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md)
and
[`distrib_d5expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d5expected_hessian.md)
when `method` is an object.

## Usage

``` r
dexpected_by_method(distrib, y, theta, scale, m, method)
```

## Arguments

- distrib:

  A univariate distribution.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned.

- scale:

  `"parameter"` or `"link"`.

- m:

  The derivative order, 1 to 5.

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md),
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md),
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)
or
[`d5expected_names()`](https://statmodels7.github.io/distributions7/reference/d5expected_names.md).

## Details

For
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
the matrix is \\E\_{ab}(\theta) = \mathbb{E}\_\theta\[f\_{ab}\]\\ and
each derivative moves the density as well as the integrand, so that the
integrand of the derivative along \\c\\ is \\\partial_c f + f\\\ell_c\\.
For
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
the matrix is \\f\_{ab}\\ evaluated at the observations, the data do not
move, and the derivative is the Leibniz derivative of \\f\_{ab}\\ alone:
\\-(\ell\_{ac}\ell_b + \ell_a\ell\_{bc})\\ at order 1 and
\\-(\ell\_{acd}\ell_b + \ell\_{ac}\ell\_{bd} + \ell\_{ad}\ell\_{bc} +
\ell_a\ell\_{bcd})\\ at order 2. The integrand \\f\_{ab}\\ is
\\\ell\_{ab}\\ for `identity = "hessian"` and \\-\ell_a\ell_b\\ for
`identity = "score"`;
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
uses the second.

The derivative of order `m` of \\f\_{ab}\\ needs observed derivatives of
order up to `m + 1` for the score identity and `m + 2` for the Hessian
identity. The package provides six orders, the sixth by the stencil of
[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md),
so the Hessian identity at `m = 5` signals an error that names
`identity = "score"`.

On the link scale and for `m` up to 4, the parameter-scale derivatives
of every order up to `m` are computed by the same method and carried
across by
[`dexpected_link_higher()`](https://statmodels7.github.io/distributions7/reference/dexpected_link_higher.md).
At `m = 5` the polynomial is evaluated with the observed derivatives on
the link scale instead, the route that
[`distrib_d5expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d5expected_hessian.md)
takes for its stencil: the identities and the moved density hold in any
coordinates, and the composition would read the parameter-scale stencil
of
[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md)
through the Faa di Bruno formula at order six.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
