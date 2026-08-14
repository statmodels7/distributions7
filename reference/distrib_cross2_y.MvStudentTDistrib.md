# Multivariate Student t Higher Mixed Response Derivatives

The three derivatives a marginal criterion reads when this family is a
prior. With \\M = \ell^{(yy)} = -c\\\Sigma^{-1} + 2d\\ww^\top\\,
\$\$\partial_a M = -c_a\Sigma^{-1} - c\\\partial_a\Sigma^{-1} + 2d_a
ww^\top + 2d\left(w_aw^\top + ww_a^\top\right),\$\$ \$\$\partial_a
\ell^{(y)} = -c_a w - c\\w_a, \qquad \partial\_{ab}\ell^{(y)} =
-c\_{ab}w - c_aw_b - c_bw_a - c\\w\_{ab},\$\$ and \\\partial\_{ab}M\\ is
the same expansion carried one order further. Every piece comes from
[`mvt_dpieces`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md).

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic.

- ...:

  Unused.

## Value

A named list, keyed by parameter or by parameter pair.

## Details

Unlike the gaussian's, this family's response Hessian depends on the
observation, so `distrib_cross2_y` and `distrib_hess_y_hess` return one
matrix per row rather than one matrix.
