# Box-Cox Transformation

Transformer for the one-parameter Box-Cox transformation \\Y =
(X^\lambda - 1)/\lambda\\ (with \\Y = \log X\\ for \\\lambda = 0\\).
Requires \\X \> 0\\.

## Usage

``` r
bc_transform(lambda)
```

## Arguments

- lambda:

  Numeric. The transformation parameter.

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
[`yj_transform`](https://statmodels7.github.io/distributions7/reference/yj_transform.md),
[`softplus_transform`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`asinh_transform`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gamma2_distrib(), bc_transform(lambda = 0.5))
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.5061537
```
