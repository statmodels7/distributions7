# Power Transformation

Transformer for \\Y = X^p\\. Inverse \\X = Y^{1/p}\\, Jacobian \\\|J\| =
\frac{1}{\|p\|}\|Y\|^{1/p - 1}\\. Fractional powers require a
non-negative support; even integer powers require a support that does
not straddle 0.

## Usage

``` r
power_transform(p = 2)
```

## Arguments

- p:

  Numeric. The exponent. Defaults to 2.

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## See also

[`log_transform`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`sqrt_transform`](https://statmodels7.github.io/distributions7/reference/sqrt_transform.md),
[`inverse_transform`](https://statmodels7.github.io/distributions7/reference/inverse_transform.md),
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
d <- transformation(gamma2_distrib(), power_transform(p = 2))
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.180447
```
