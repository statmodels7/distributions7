# Exponential Transformation

Transformer for \\Y = e^X\\: maps the real line to \\(0, \infty)\\.
Inverse \\X = \log(Y)\\, Jacobian \\\|J\| = 1/Y\\.

## Usage

``` r
exp_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#> [1] 0.3989423
```
