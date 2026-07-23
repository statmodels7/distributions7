# Inverse Hyperbolic Sine Transformation

Transformer for \\Y = \text{asinh}(X)\\, a log-like transformation that
handles zero and negative values. Inverse \\X = \sinh(Y)\\, Jacobian
\\\|J\| = \cosh(Y)\\.

## Usage

``` r
asinh_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.
