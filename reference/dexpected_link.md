# The Derivatives of the Link-Scale Expected Information

Leibniz's rule on \\F\_{ab} = E\_{ab}\\h_a' h_b'\\, written once for
every family.

## Usage

``` r
dexpected_link(params, E, d1, d2, h, order)
```

## Arguments

- params:

  The parameter names.

- E:

  The expected information on the parameter scale.

- d1:

  Its first derivatives, keyed as
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md).

- d2:

  Its second derivatives, keyed as
  [`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md),
  or `NULL` at order 1.

- h:

  The inverse link's derivatives, from
  [`inverse_link_derivs()`](https://statmodels7.github.io/distributions7/reference/inverse_link_derivs.md),
  to order `order + 1`.

- order:

  `1L` or `2L`.

## Value

A named list on the link scale.

## Details

Each parameter moves with its own coordinate only, so with \\u\_{ab} =
h_a' h_b'\\, \$\$\partial_c u = \[a{=}c\]\\h_a'' h_b' + \[b{=}c\]\\h_a'
h_b'',\$\$ \$\$\partial\_{cd} u = \[a{=}c{=}d\]\\h_a''' h_b' +
(\[a{=}c\]\[b{=}d\] + \[b{=}c\]\[a{=}d\])\\h_a'' h_b'' +
\[b{=}c{=}d\]\\h_a' h_b''',\$\$ and \$\$\partial_c F = (\partial_c
E)\\h_c' u + E\\\partial_c u,\$\$ \$\$\partial\_{cd} F = (\partial\_{cd}
E)\\h_c' h_d' u + \[c{=}d\]\\(\partial_c E)\\h_c'' u + (\partial_c
E)\\h_c'\\\partial_d u + (\partial_d E)\\h_d'\\\partial_c u +
E\\\partial\_{cd} u.\$\$
