# Logit Transformation

Transformer for \\Y = \text{logit}(X)\\: maps \\(0, 1)\\ to the real
line. Inverse \\X = \text{plogis}(Y)\\, Jacobian \\\|J\| =
\text{dlogis}(Y)\\.

## Usage

``` r
logit_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(beta_distrib(), logit_transform())
distrib_pdf(d, 0, list(mu = 0.4, phi = 5))
#> [1] 0.375
```
