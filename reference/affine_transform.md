# Affine (Location-Scale) Transformation

Transformer for \\Y = \text{loc} + \text{scale} \cdot X\\. Inverse \\X =
(Y - \text{loc})/\text{scale}\\, Jacobian \\\|J\| =
1/\|\text{scale}\|\\.

## Usage

``` r
affine_transform(loc = 0, scale = 1)
```

## Arguments

- loc:

  Numeric. The location shift. Defaults to 0.

- scale:

  Numeric. The scale multiplier (non-zero). Defaults to 1.

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
[`asinh_transform()`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform()`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform()`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gaussian1_distrib(), affine_transform(loc = 1, scale = 2))
distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#> [1] 0.1994711
```
