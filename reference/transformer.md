# S7 Class for Variable Transformers

A `transformer` object defines the mathematical rules for transforming a
random variable \\Y = g(X)\\. It is used as input to
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md).

## Usage

``` r
transformer(
  name = character(0),
  trans_fun = function() NULL,
  trans_inv = function() NULL,
  trans_abs_jac = function() NULL,
  trans_inv_hessian = function() NULL,
  grad_log_jac = function() NULL,
  hess_log_jac = function() NULL,
  bounds_fun = function() NULL,
  valid_support = function() NULL,
  decreasing = logical(0)
)
```

## Arguments

- name:

  A character string identifying the transformation.

- trans_fun:

  The forward transformation \\y = g(x)\\.

- trans_inv:

  The inverse transformation \\x = g^{-1}(y)\\.

- trans_abs_jac:

  The absolute Jacobian of the inverse transformation \\\|J(y)\| =
  \|dx/dy\|\\; must accept a `log` argument.

- trans_inv_hessian:

  The second derivative of the inverse transformation \\d^2x/dy^2\\.

- grad_log_jac:

  The first derivative of \\\log\|J(y)\|\\ with respect to \\y\\.

- hess_log_jac:

  The second derivative of \\\log\|J(y)\|\\ with respect to \\y\\.

- bounds_fun:

  Maps the original support bounds to the transformed ones.

- valid_support:

  Checks whether a support is compatible with the transformation.

- decreasing:

  Logical; `TRUE` for monotonically decreasing transformations.

## Value

An object of class `transformer`.

## Methods

No method dispatches on this class: a `transformer` is a description of
a change of variables, not a distribution. It is consumed by
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which returns a
[`TransformedDistrib`](https://statmodels7.github.io/distributions7/reference/TransformedDistrib.md)
carrying the full set of distribution methods.

Ready-made transformers:
[`log_transform`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`sqrt_transform`](https://statmodels7.github.io/distributions7/reference/sqrt_transform.md),
[`inverse_transform`](https://statmodels7.github.io/distributions7/reference/inverse_transform.md),
[`power_transform`](https://statmodels7.github.io/distributions7/reference/power_transform.md),
[`bc_transform`](https://statmodels7.github.io/distributions7/reference/bc_transform.md),
[`yj_transform`](https://statmodels7.github.io/distributions7/reference/yj_transform.md),
[`softplus_transform`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`asinh_transform`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md).

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
tr <- log_transform()
S7::S7_inherits(tr, transformer)
#> [1] TRUE
tr@name
#> [1] "log"

# a transformer is consumed by transformation(), which is where it acts
distrib_pdf(transformation(gamma2_distrib(), tr), 0, list(mu = 2, sigma2 = 1))
#> [1] 0.3608941
```
