# The Terms of a Derivative of the Expectation of a Product of Two Scores

Expands \\\partial_S E\[\ell_a \ell_b\]\\ as the expectation of a
polynomial in the derivatives of the log-density \\\ell\\.

## Usage

``` r
d34e_plan(P, a, b, S)
```

## Arguments

- P:

  The parameter names.

- a, b:

  The indices of the pair.

- S:

  The integer indices of the derivative, with repetition.

## Value

A list with `coef`, the coefficient of each term, `factors`, a list of
the keys of each term's factors, and `blocks`, a named list giving the
parameter names of each key.

## Details

With \\f = e^{\ell}\\ the derivative of the density in the parameters of
a multiset \\R\\ is \\f\\ times the complete Bell polynomial \\Y_R\\ of
the derivatives of \\\ell\\, the sum over the set partitions of \\R\\ of
the products of \\\ell^{(B)}\\ over the blocks. Then \$\$\partial_S
E\[\ell_a\ell_b\] = E\Big\[\sum\_{T \subseteq S}
\partial_T(\ell_a\ell_b)\\ Y\_{S \setminus T}\Big\],\qquad
\partial_T(\ell_a\ell_b) = \sum\_{T = T_1 \sqcup T_2}
\ell\_{aT_1}\\\ell\_{bT_2},\$\$ where the sums run over the subsets of
the positions of \\S\\. The terms are collected by their sorted factors.
The highest derivative of \\\ell\\ that appears is of order \\\|S\| +
1\\.
