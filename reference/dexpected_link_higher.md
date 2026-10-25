# The Derivatives of the Link-Scale Expected Information, Any Order

Applies Leibniz's rule and the diagonal Faa di Bruno formula to
\\F\_{ab}(\eta) = E\_{ab}(\theta(\eta))\\h_a'(\eta_a)\\h_b'(\eta_b)\\,
at orders one to four.

## Usage

``` r
dexpected_link_higher(params, nat, h, order)
```

## Arguments

- params:

  The parameter names.

- nat:

  A list: `nat[[1]]` the expected information on the parameter scale,
  keyed as
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
  and `nat[[m + 1]]` its derivatives of order `m`, keyed as
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
  [`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md),
  [`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
  and
  [`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

- h:

  The inverse link's derivatives, from
  [`inverse_link_derivs()`](https://statmodels7.github.io/distributions7/reference/inverse_link_derivs.md),
  to order `order + 1`.

- order:

  `1L` to `4L`.

## Value

A named list on the link scale, keyed as the order's names.

## Details

For a multiset \\S\\ of coordinates, with \\G\_{ab} =
E\_{ab}\circ\theta\\ and \\u\_{ab} = h_a' h_b'\\, \$\$\partial_S F\_{ab}
= \sum\_{T \subseteq S} \binom{S}{T}\\ \partial_T G\_{ab}\\\partial\_{S
\setminus T} u\_{ab},\$\$ where the sum runs over the sub-multisets of
\\S\\ and the coefficient is the product over the distinct indices of
the binomial coefficients of their multiplicities. Each parameter moves
with its own coordinate only, so \$\$\partial_T G\_{ab} = \sum\_{j}
\partial^{J} E\_{ab} \prod_i B\_{m_i, j_i}(h_i', h_i'', \dots),\$\$ with
\\m_i\\ the multiplicity of index \\i\\ in \\T\\, \\1 \le j_i \le m_i\\,
\\J\\ the multiset holding each \\i\\ \\j_i\\ times, and \\B\\ the
partial Bell polynomials of
[`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md).
The factor \\\partial_R u\_{ab}\\ vanishes unless \\R\\ names only \\a\\
and \\b\\; then it is \\h_a^{(1+r_a)} h_b^{(1+r_b)}\\ for \\a \ne b\\
and \\\sum_k \binom{r}{k} h_a^{(1+k)} h_a^{(1+r-k)}\\ for \\a = b\\.
