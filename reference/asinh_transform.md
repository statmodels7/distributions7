# Inverse Hyperbolic Sine Transformation

Transformer for \\Y = \text{asinh}(X)\\, a log-like transformation that
handles zero and negative values. Inverse \\X = \sinh(Y)\\, Jacobian
\\\|J\| = \cosh(Y)\\.

## Usage

``` r
asinh_transform()
```

## Value

A
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## See also

[`log_transform()`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform()`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`sqrt_transform()`](https://statmodels7.github.io/distributions7/reference/sqrt_transform.md),
[`inverse_transform()`](https://statmodels7.github.io/distributions7/reference/inverse_transform.md),
[`power_transform()`](https://statmodels7.github.io/distributions7/reference/power_transform.md),
[`bc_transform()`](https://statmodels7.github.io/distributions7/reference/bc_transform.md),
[`yj_transform()`](https://statmodels7.github.io/distributions7/reference/yj_transform.md),
[`softplus_transform()`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`logit_transform()`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform()`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform()`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gaussian1_distrib(), asinh_transform())
distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#> [1] 0.3086008
```
