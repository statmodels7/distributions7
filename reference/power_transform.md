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
