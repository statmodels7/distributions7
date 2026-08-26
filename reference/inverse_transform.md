# Reciprocal (Inverse) Transformation

Transformer for \\Y = 1/X\\ (monotonically decreasing on a support not
containing 0). Inverse \\X = 1/Y\\, Jacobian \\\|J\| = 1/Y^2\\.

## Usage

``` r
inverse_transform()
```

## Value

A
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## See also

[`log_transform()`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform()`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`sqrt_transform()`](https://statmodels7.github.io/distributions7/reference/sqrt_transform.md),
[`power_transform()`](https://statmodels7.github.io/distributions7/reference/power_transform.md),
[`bc_transform()`](https://statmodels7.github.io/distributions7/reference/bc_transform.md),
[`yj_transform()`](https://statmodels7.github.io/distributions7/reference/yj_transform.md),
[`softplus_transform()`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`asinh_transform()`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform()`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform()`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform()`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gamma2_distrib(), inverse_transform())
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.3608941
```
