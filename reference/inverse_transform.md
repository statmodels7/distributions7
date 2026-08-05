# Reciprocal (Inverse) Transformation

Transformer for \\Y = 1/X\\ (monotonically decreasing on a support not
containing 0). Inverse \\X = 1/Y\\, Jacobian \\\|J\| = 1/Y^2\\.

## Usage

``` r
inverse_transform()
```

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(gamma2_distrib(), inverse_transform())
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.3608941
```
