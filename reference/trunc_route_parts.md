# The Retained Mass of a Truncated Family and the Sums of Its Information

Computes, for a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family whose parent has a scalar route, the retained mass \\Z\\, its
derivatives, and the sums that the truncated expected information and
its first derivatives read, by the same arithmetic as the compiled
registry.

## Usage

``` r
trunc_route_parts(distrib, y, theta, what, rows = NULL)
```

## Arguments

- distrib:

  A truncated family.

- y:

  The observations, read for their number.

- theta:

  A named list of parameter values.

- what:

  One of `"z"`, `"grad"`, `"hess"`, `"info"` and `"dinfo"`, each adding
  quantities to those of the previous one.

- rows:

  `NULL`, or the observations to compute; the function sets it when it
  computes each distinct parameter vector once.

## Value

`NULL` when the family has no route; otherwise a list with `Z`, `Zi`
(named by parameter), and, as `what` requires, `Zij` and `S` (named as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md))
and `D` (named as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)),
each a vector over the observations.

## Details

With \\f\\, \\s_i\\ and \\\ell\_{ij}\\ the parent's density, score and
Hessian, the quantities are \$\$Z = \int f, \quad Z_i = \int f s_i,
\quad Z\_{ij} = \int f (\ell\_{ij} + s_i s_j), \quad S\_{ij} = \int f
s_i s_j,\$\$ \$\$D\_{ijc} = \int f (s_i s_j s_c + \ell\_{ic} s_j + s_i
\ell\_{jc}),\$\$ over the retained set, so that \\E_T\[s_i s_j\] =
S\_{ij}/Z\\ and \\D\_{ijc}\\ is the derivative of \\\int f s_i s_j\\ in
parameter \\c\\.

For a discrete parent the integrals are sums over support points.
`trunc_rule_cpp()` returns, for each observation, the points they run
over: the retained points when they are finitely many, and otherwise a
series over the retained points stopped by a rule on the masses alone.
For the density, the score and the Hessian on an unbounded support, the
points removed below `lower` are used instead when their mass is at most
one half, and the retained sums are \\1 - \sum f\\, \\-\sum f s_i\\ and
\\-\sum f(\ell\_{ij} + s_i s_j)\\.

For a continuous parent, \\Z\\, \\Z_i\\ and \\Z\_{ij}\\ are the
differences of the compiled distribution function and its derivatives at
the two truncation points (`trunc_cont_ends_cpp()`), \\Z\\ in the
survival function when the lower point lies above the median; \\S\\ and
\\D\\ are taken by the rule of `trunc_cont_rule_cpp()`, whose nodes
depend on the interval and on the parent's center, scale and kinks only.

The sums accumulate in long double in the registry's order, so that each
diagonal quantity is the registry's to the bit.
