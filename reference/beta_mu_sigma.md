# The Mean Proportion and Dispersion of a Beta in Its Shapes

Returns \\\mu = \alpha/(\alpha + \beta)\\ as \\1/(1 + \beta/\alpha)\\,
so that it is finite where the sum of the shapes overflows, and \\\sigma
= 1/(\alpha + \beta)\\.

## Usage

``` r
beta_mu_sigma(a, b)
```

## Arguments

- a, b:

  The two shapes, positive numeric vectors.

## Value

A list with components `mu` and `sigma`.

## See also

[`betabinom_moments()`](https://statmodels7.github.io/distributions7/reference/betabinom_moments.md),
which takes them.
