# The Key of One Component of the Expected Information's Derivative

The name under which
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
returns \\\partial\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\.

## Usage

``` r
dexpected_key(params, a, b, k)
```

## Arguments

- params:

  A character vector of parameter names, in the family's order.

- a, b:

  Indices into `params`; their order does not matter, the component
  being symmetric in them.

- k:

  The index of the parameter differentiated in, which does matter.

## Value

A single string.

## See also

[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)

## Examples

``` r
dexpected_key(c("mu", "sigma"), 1, 2, 2)
#> [1] "mu_sigma_sigma"
```
