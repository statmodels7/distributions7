# S7 Class for Variable Transformers

A `transformer` records the mathematical rules of a change of variables
\\Y = g(X)\\: the map, its inverse, the Jacobian of that inverse and the
two derivatives of the Jacobian's logarithm, together with the
bookkeeping that says which supports the map applies to and what it does
to them. It carries no distribution and answers no distribution generic;
it is consumed by
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which returns a
[TransformedDistrib](https://statmodels7.github.io/distributions7/reference/TransformedDistrib.md).

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

  A single string identifying the transformation. It composes the
  result's `distrib_name` as `"name(parent)"` unless
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
  is given a `new_name`.

- trans_fun:

  The forward map \\y = g(x)\\, a function of one numeric vector
  returning one of the same length.

- trans_inv:

  The inverse map \\x = g^{-1}(y)\\, likewise.

- trans_abs_jac:

  The absolute Jacobian of the inverse, \\\lvert J(y)\rvert = \lvert
  dx/dy\rvert\\. It MUST accept a `log` argument and return the
  logarithm when it is `TRUE`.

- trans_inv_hessian:

  The second derivative of the inverse map, \\d^2x/dy^2\\.

- grad_log_jac:

  The first derivative of \\\log\lvert J(y)\rvert\\ with respect to
  \\y\\.

- hess_log_jac:

  The second derivative of \\\log\lvert J(y)\rvert\\ with respect to
  \\y\\.

- bounds_fun:

  A function of the parent's `bounds` returning the transformed support,
  a numeric vector of length 2.

- valid_support:

  A function of the parent's `bounds` returning a single logical:
  whether the map applies there. It is what
  `transformation(gaussian1_distrib(), log_transform())` fails on.

- decreasing:

  Logical of length 1. `TRUE` for a monotonically decreasing map, which
  swaps the tails in the distribution and quantile functions.

## Value

An S7 object of class `transformer`, carrying the ten properties above
and nothing else.

## What a map must be to be one

\\g\\ has to be BIJECTIVE on the parent's support. A map that is two to
one, such as the absolute value, has no inverse to carry a density
through and cannot be a transformer at all;
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
handles that case by adding the two preimages instead.

## Writing your own

The ten properties are all required, and `trans_abs_jac` must accept a
`log` argument, because
[`distrib_pdf.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.TransformedDistrib.md)
works on the log scale throughout. `valid_support` takes the parent's
`bounds` and answers whether the map applies; `bounds_fun` takes the
same and returns the transformed support. Nothing validates the
derivatives, so
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
on the resulting distribution is what catches a transcription error in
them.

## Notation

\\g\\ is the transformation, \\J\\ the Jacobian of its inverse, \\X\\
the parent's variable and \\Y = g(X)\\ the transformed one.

## Methods

None. A `transformer` is a description of a change of variables, not a
distribution. It is consumed by
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which returns a
[TransformedDistrib](https://statmodels7.github.io/distributions7/reference/TransformedDistrib.md)
carrying the full set of distribution methods.

The twelve ready-made transformers:
[`log_transform()`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform()`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`sqrt_transform()`](https://statmodels7.github.io/distributions7/reference/sqrt_transform.md),
[`inverse_transform()`](https://statmodels7.github.io/distributions7/reference/inverse_transform.md),
[`power_transform()`](https://statmodels7.github.io/distributions7/reference/power_transform.md),
[`bc_transform()`](https://statmodels7.github.io/distributions7/reference/bc_transform.md),
[`yj_transform()`](https://statmodels7.github.io/distributions7/reference/yj_transform.md),
[`softplus_transform()`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`asinh_transform()`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform()`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform()`](https://statmodels7.github.io/distributions7/reference/expit_transform.md)
and
[`affine_transform()`](https://statmodels7.github.io/distributions7/reference/affine_transform.md).

## See also

[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which consumes one,
[TransformedDistrib](https://statmodels7.github.io/distributions7/reference/TransformedDistrib.md)
for what it produces, and
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
for a map that cannot be a transformer.

## Examples

``` r
tr <- log_transform()
S7::S7_inherits(tr, transformer)
#> [1] TRUE
tr@name
#> [1] "log"
S7::prop_names(tr)
#>  [1] "name"              "trans_fun"         "trans_inv"        
#>  [4] "trans_abs_jac"     "trans_inv_hessian" "grad_log_jac"     
#>  [7] "hess_log_jac"      "bounds_fun"        "valid_support"    
#> [10] "decreasing"       

# The map, its inverse and the Jacobian of that inverse.
c(forward = tr@trans_fun(exp(1)), inverse = tr@trans_inv(1),
  log_jacobian = tr@trans_abs_jac(1, log = TRUE))
#>      forward      inverse log_jacobian 
#>     1.000000     2.718282     1.000000 

# It says which supports it applies to, and what it does to them.
c(on_positives = tr@valid_support(c(0, Inf)),
  on_the_line = tr@valid_support(c(-Inf, Inf)))
#> on_positives  on_the_line 
#>         TRUE        FALSE 
tr@bounds_fun(c(0, Inf))
#> [1] -Inf  Inf

# And a decreasing map declares itself, which swaps the tails downstream.
vapply(list(log_transform(), inverse_transform(), expit_transform()),
       function(z) z@decreasing, TRUE)
#> [1] FALSE  TRUE FALSE

# A transformer acts only through transformation().
distrib_pdf(transformation(gamma2_distrib(), tr), 0,
            list(mu = 2, sigma2 = 1))
#> [1] 0.3608941
```
