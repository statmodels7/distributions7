# Softplus Transformation

Returns the transformer for the inverse-softplus map \\Y =
\frac{1}{a}\log(e^{aX} - 1)\\, which sends \\(0, \infty)\\ to the real
line and is evaluated as \\X + \frac{1}{a}\log(1 - e^{-aX})\\ so that it
does not overflow at a large argument. The inverse (softplus) is \\X =
\frac{1}{a}\log(1 + e^{aY})\\ and the Jacobian is \\\|J\| =
\text{plogis}(aY)\\.

## Usage

``` r
softplus_transform(a = 1)
```

## Arguments

- a:

  The numeric scale parameter, which must be positive. Defaults to 1.

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
[`asinh_transform()`](https://statmodels7.github.io/distributions7/reference/asinh_transform.md),
[`logit_transform()`](https://statmodels7.github.io/distributions7/reference/logit_transform.md),
[`expit_transform()`](https://statmodels7.github.io/distributions7/reference/expit_transform.md),
[`affine_transform()`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- transformation(gamma2_distrib(), softplus_transform(a = 1))
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.3193671
```
