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

## Examples

``` r
d <- transformation(gamma_distrib(), sqrt_transform())
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.7217882
```
