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
