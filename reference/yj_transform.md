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
