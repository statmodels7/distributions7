# The Expected Information's Derivatives Through a Reparametrization

Carries the parent's \\\partial E/\partial\theta\\ and \\\partial^2
E/\partial\theta^2\\ into the new coordinates \\\phi\\.

## Usage

``` r
reparam_dexpected(distrib, y, theta, order, threads = 1L)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  The response.

- theta:

  A named list of the new parameters.

- order:

  `1L` or `2L`.

- threads:

  The thread count passed to the parent.

## Value

A named list on the parameter scale, keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

The expected information transforms as a tensor, with no term in the
map's second derivative because the score has mean zero: \\\tilde
E\_{ab} = \sum\_{ij} E\_{ij} J\_{ia} J\_{jb}\\, with \\J\_{ia} =
\partial\theta_i/\partial\phi_a\\. Writing \\J\_{ia,c}\\ and
\\J\_{ia,cd}\\ for the map's second and third derivatives,
\$\$\partial_c\tilde E\_{ab} = \sum\_{ij}\Big\[\sum_k \partial_k E\_{ij}
J\_{kc}\Big\] J\_{ia}J\_{jb} + \sum\_{ij} E\_{ij}(J\_{ia,c}J\_{jb} +
J\_{ia}J\_{jb,c}),\$\$ and \\\partial\_{cd}\tilde E\_{ab}\\ follows by
Leibniz's rule once more, reading the parent's second derivative and the
map's third. The map's partials are the keyed tables of
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md),
hand-written for the shipped reparametrizations, so every term is
analytic wherever the parent's derivatives are.

## See also

[`reparam_chain()`](https://statmodels7.github.io/distributions7/reference/reparam_chain.md),
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md)
