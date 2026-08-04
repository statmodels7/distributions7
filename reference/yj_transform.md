# Yeo-Johnson Transformation

Transformer for the Yeo-Johnson transformation, extending Box-Cox to
negative values: Box-Cox of \\X+1\\ (parameter \\\lambda\\) for \\X \ge
0\\, and negated Box-Cox of \\\|X\|+1\\ (parameter \\2-\lambda\\) for
\\X \< 0\\.

## Usage

``` r
yj_transform(lambda)
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
d <- transformation(gaussian_distrib(), yj_transform(lambda = 0.5))
distrib_pdf(d, 1, list(mu = 0, sigma = 1))
#> [1] 0.2739736
```
