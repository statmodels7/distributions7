# Logit Transformation

Transformer for \\Y = \text{logit}(X)\\: maps \\(0, 1)\\ to the real
line. Inverse \\X = \text{plogis}(Y)\\, Jacobian \\\|J\| =
\text{dlogis}(Y)\\.

## Usage

``` r
logit_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## See also

[`log_transform`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`sqrt_transform`](https://statmodels7.github.io/distributions7/reference/sqrt_transform.md),
[`inverse_transform`](https://statmodels7.github.io/distributions7/reference/inverse_transform.md),
[`power_transform`](https://statmodels7.github.io/distributions7/reference/power_transform.md),
[`bc_transform`](https://statmodels7.github.io/distributions7/reference/bc_transform.md),
[`yj_transform`](https://statmodels7.github.io/distributions7/reference/yj_transform.md),
[`softplus_transform`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`asinh_transform`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`expit_transform`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(beta1_distrib(), logit_transform())
distrib_pdf(d, 0, list(mu = 0.4, phi = 5))
#> [1] 0.375
```
