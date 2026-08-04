# Logarithmic Transformation

Transformer for \\Y = \log(X)\\: maps \\(0, \infty)\\ to the real line.
Inverse \\X = e^Y\\, Jacobian \\\|J\| = e^Y\\.

## Usage

``` r
log_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(gamma_distrib(), log_transform())
distrib_pdf(d, 0, list(mu = 2, sigma2 = 1))
#> [1] 0.3608941
```
