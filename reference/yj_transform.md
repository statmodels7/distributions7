# Yeo-Johnson Transformation

Transformer for the Yeo-Johnson transformation, extending Box-Cox to
negative values: Box-Cox of \\X+1\\ (parameter \\\lambda\\) for \\X \ge
0\\, and negated Box-Cox of \\\|X\|+1\\ (parameter \\2-\lambda\\) for
\\X \< 0\\.

## Usage

``` r
yj_transform(lambda)
```

## Arguments

- lambda:

  Numeric. The transformation parameter.

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
[`softplus_transform()`](https://statmodels7.github.io/distributions7/reference/softplus_transform.md),
[`asinh_transform()`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform()`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform()`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform()`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gaussian1_distrib(), yj_transform(lambda = 0.5))
distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#> [1] 0.2739736
```
