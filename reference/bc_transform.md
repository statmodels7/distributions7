# Box-Cox Transformation

Transformer for the one-parameter Box-Cox transformation \\Y =
(X^\lambda - 1)/\lambda\\ (with \\Y = \log X\\ for \\\lambda = 0\\).
Requires \\X \> 0\\.

## Usage

``` r
bc_transform(lambda)
```

## Arguments

- lambda:

  Numeric. The transformation parameter.

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(gamma_distrib(), bc_transform(lambda = 0.5))
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.5061537
```
