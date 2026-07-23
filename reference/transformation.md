# Apply a Variable Transformation to a Distribution Object

Creates a new distribution object for \\Y = g(X)\\, where \\X\\ follows
an existing **continuous** distribution and \\g\\ is a bijective
transformation described by a
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
object.

## Usage

``` r
transformation(distrib, transformer)
```

## Arguments

- distrib:

  An object inheriting from `continuous_distrib`.

- transformer:

  A
  [`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md)
  object (e.g.
  [`log_transform()`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
  [`bc_transform`](https://statmodels7.github.io/distributions7/reference/bc_transform.md)`(0.5)`,
  [`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md)`(1, 2)`).

## Value

An S7 object of class `TransformedDistrib` (inheriting from
`continuous_distrib`).

## Details

The density follows the change-of-variables formula \$\$f_Y(y) =
f_X(g^{-1}(y)) \cdot \left\|\dfrac{d}{dy} g^{-1}(y)\right\|\$\$ CDF,
quantiles and RNG are obtained by mapping through \\g\\ (with tails
swapped for decreasing transformations). Since \\g\\ does not depend on
the parameters, the score, observed Hessian and expected Hessian
coincide with the parent's, evaluated at \\x = g^{-1}(y)\\. Moments are
available numerically via
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).

## See also

[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md),
[`log_transform`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
[`exp_transform`](https://statmodels7.github.io/distributions7/reference/exp_transform.md),
[`affine_transform`](https://statmodels7.github.io/distributions7/reference/affine_transform.md),
[`bc_transform`](https://statmodels7.github.io/distributions7/reference/bc_transform.md),
[`yj_transform`](https://statmodels7.github.io/distributions7/reference/yj_transform.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# A lognormal built by transformation, equal to lognormal_distrib()
logn <- transformation(gaussian_distrib(), exp_transform())
distrib_pdf(logn, 2, list(mu = 0, sigma = 1))
dlnorm(2, 0, 1)
} # }
```
