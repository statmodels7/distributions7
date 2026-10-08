# The Key of One Component of the Expected Information's Third Derivative

The name under which
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md)
returns
\\\partial^3\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\\partial\theta_d\\
\partial\theta_e\\.

## Usage

``` r
d3expected_key(params, a, b, c, d, e)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

- a, b:

  Indices of the information's pair; their order does not matter.

- c, d, e:

  Indices of the parameters differentiated in; their order does not
  matter either.

## Value

A single string.

## See also

[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)

## Examples

``` r
d3expected_key(c("mu", "sigma"), 1, 1, 2, 1, 2)
#> [1] "mu_mu_mu_sigma_sigma"
```
