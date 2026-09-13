# The Key of One Component of the Expected Information's Second Derivative

The name under which
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
returns
\\\partial^2\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\\partial\theta_d\\.

## Usage

``` r
d2expected_key(params, a, b, c, d)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

- a, b:

  Indices of the information's pair; their order does not matter.

- c, d:

  Indices of the parameters differentiated in; their order does not
  matter either.

## Value

A single string.

## See also

[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md),
[`dexpected_key()`](https://statmodels7.github.io/distributions7/reference/dexpected_key.md)

## Examples

``` r
d2expected_key(c("mu", "sigma"), 1, 1, 2, 2)
#> [1] "mu_mu_sigma_sigma"
```
