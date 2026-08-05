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

## Examples

``` r
d <- transformation(gaussian1_distrib(), asinh_transform())
distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#> [1] 0.3086008
```
