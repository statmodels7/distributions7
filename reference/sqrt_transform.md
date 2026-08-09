# Square Root Transformation

Transformer for \\Y = \sqrt{X}\\ (requires \\X \ge 0\\). Inverse \\X =
Y^2\\, Jacobian \\\|J\| = 2Y\\.

## Usage

``` r
sqrt_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## See also

[`log_transform`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`inverse_transform`](https://statmodels7.github.io/distributions7/reference/inverse_transform.md),
[`power_transform`](https://statmodels7.github.io/distributions7/reference/power_transform.md),
[`bc_transform`](https://statmodels7.github.io/distributions7/reference/bc_transform.md),
[`yj_transform`](https://statmodels7.github.io/distributions7/reference/yj_transform.md),
[`softplus_transform`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`asinh_transform`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gamma2_distrib(), sqrt_transform())
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.7217882
```
