# Expit (Sigmoid) Transformation

Transformer for \\Y = \text{plogis}(X)\\: maps the real line to \\(0,
1)\\. Inverse \\X = \text{logit}(Y)\\, Jacobian \\\|J\| = 1/(Y(1-Y))\\.

## Usage

``` r
expit_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(gaussian_distrib(), expit_transform())
distrib_pdf(d, 0.5, list(mu = 0, sigma = 1))
#> [1] 1.595769
```
