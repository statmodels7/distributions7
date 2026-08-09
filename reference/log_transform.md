# Logarithmic Transformation

Transformer for \\Y = \log(X)\\: maps \\(0, \infty)\\ to the real line.
Inverse \\X = e^Y\\, Jacobian \\\|J\| = e^Y\\.

## Usage

``` r
log_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## See also

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
[`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gamma2_distrib(), log_transform())
distrib_pdf(d, 0, list(mu = 2, sigma2 = 1))
#> [1] 0.3608941
```
