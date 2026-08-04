# Softplus Transformation

Transformer for the inverse-softplus map \\Y = \frac{1}{a}\log(e^{aX} -
1)\\, sending \\(0, \infty)\\ to the real line. Inverse (softplus) \\X =
\frac{1}{a}\log(1 + e^{aY})\\, Jacobian \\\|J\| = \text{plogis}(aY)\\.

## Usage

``` r
softplus_transform(a = 1)
```

## Arguments

- a:

  Numeric. Positive scale parameter. Defaults to 1.

## Value

A
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Examples

``` r
d <- transformation(gamma_distrib(), softplus_transform(a = 1))
distrib_pdf(d, 1, list(mu = 2, sigma2 = 1))
#> [1] 0.3193671
```
