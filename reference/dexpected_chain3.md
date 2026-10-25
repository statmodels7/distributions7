# The Third Derivative of the Expected Information Through a Reparametrization

Computes the third derivative of the expected information in the new
coordinates \\\phi\\ from the parent's expected information and its
first three derivatives, and the map's partial derivatives up to the
fourth order.

## Usage

``` r
dexpected_chain3(Pp, Pn, E, d1, d2, d3, maps, n)
```

## Arguments

- Pp, Pn:

  The parent's parameter names and the new ones.

- E, d1, d2, d3:

  Named lists keyed as
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
  [`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
  and
  [`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
  on the parent's parameters.

- maps:

  One keyed table per parent parameter, as
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  returns them, to the order 4.

- n:

  The number of observations, to which every component is recycled.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
on the new parameters.

## Details

With \\J\_{ia} = \partial\theta_i/\partial\phi_a\\, the expected
information transforms as \\\tilde E\_{ab} = \sum\_{ij} G\_{ij} J\_{ia}
J\_{jb}\\, where \\G\_{ij}(\phi) = E\_{ij}(\theta(\phi))\\, and there is
no term in the map's second derivative because the score has mean zero.
For a multiset \\S = \\c, d, e\\\\ of new coordinates, Leibniz's rule
over the three factors gives \$\$\partial_S \tilde E\_{ab} = \sum\_{ij}
\sum\_{S = S_G \cup S_A \cup S_B} \partial\_{S_G} G\_{ij}\\
\theta\_{i,\\a S_A}\\ \theta\_{j,\\b S_B},\$\$ the sum running over the
\\3^3\\ assignments of the three positions to the factors, with
\\\theta\_{i,\\a S_A}\\ the partial derivative of \\\theta_i\\ in
\\\phi_a\\ and the coordinates of \\S_A\\. The derivatives of
\\G\_{ij}\\ are those of a composition of the parent's \\E\_{ij}\\ with
the map, written out to third order: \$\$\partial_c G = \sum_k
E\_{,k}\\\theta\_{k,c},\$\$ \$\$\partial\_{cd} G = \sum\_{kl}
E\_{,kl}\\\theta\_{k,c}\theta\_{l,d} + \sum_k
E\_{,k}\\\theta\_{k,cd},\$\$ \$\$\partial\_{cde} G = \sum\_{klm}
E\_{,klm}\\\theta\_{k,c}\theta\_{l,d} \theta\_{m,e} + \sum\_{kl}
E\_{,kl}\big(\theta\_{k,cd}\theta\_{l,e} + \theta\_{k,ce}\theta\_{l,d} +
\theta\_{k,de}\theta\_{l,c}\big) + \sum_k E\_{,k}\\\theta\_{k,cde},\$\$
where \\E\_{,k}\\, \\E\_{,kl}\\ and \\E\_{,klm}\\ are the parent's
first, second and third derivatives of \\E\_{ij}\\. The map is read up
to the fourth order, \\\theta\_{i,\\a S_A}\\ with \\\|S_A\| = 3\\ being
the highest partial.

## See also

[`dexpected_chain()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain.md)
for the orders one and two.
