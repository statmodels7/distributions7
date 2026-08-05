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

## Examples

``` r
d <- transformation(gamma2_distrib(), power_transform(p = 2))
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.180447
```
