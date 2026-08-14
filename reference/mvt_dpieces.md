# The Pieces Every Mixed Response Derivative of the t Is Written In

The response derivatives of this family are \\\ell^{(y)} = -c\\w\\ and
\\\ell^{(yy)} = -c\\\Sigma^{-1} + 2d\\ww^\top\\ with \\c = (\nu+p)/s\\,
\\d = (\nu+p)/s^2\\ and \\s = \nu + q\\, so every derivative in the
parameters is assembled from the first and second derivatives of four
things: \\s\\, \\w\\, \\\Sigma^{-1}\\ and the scalars \\c\\ and \\d\\
that follow from \\s\\.

## Usage

``` r
mvt_dpieces(distrib, y, theta)
```

## Arguments

- distrib:

  A
  [`MvStudentTDistrib`](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

## Value

A list of the quantities above, indexed by parameter POSITION – the
means, then the matrix parameter's free values, then \\\nu\\ – with
`pair()` a function returning the second-order pieces of a pair.

## Details

Writing \\A_k\\ for the matrix parameter's first derivatives, \\u_j\\
for the \\j\\th column of \\\Sigma^{-1}\\ and \\P_k\\ for
\\\Sigma^{-1}A_k\Sigma^{-1}\\, the first derivatives are
\$\$\partial\_{\mu_j} s = -2w_j, \quad \partial\_{\eta_k} s = -w^\top
A_kw, \quad \partial\_\nu s = 1,\$\$ \$\$\partial\_{\mu_j} w = -u_j,
\quad \partial\_{\eta_k} w = -\Sigma^{-1}A_kw, \quad \partial\_\nu w =
0, \qquad \partial\_{\eta_k}\Sigma^{-1} = -P_k,\$\$ and the second ones
vanish except for \$\$\partial\_{\mu_j\mu_l} s = 2(\Sigma^{-1})\_{jl},
\quad \partial\_{\mu_j\eta_k} s = 2(\Sigma^{-1}A_kw)\_j, \quad
\partial\_{\mu_j\eta_k} w = P_ke_j,\$\$ and, with \\B\_{kl} =
A_k\Sigma^{-1}A_l + A_l\Sigma^{-1}A_k - A\_{kl}\\,
\$\$\partial\_{\eta_k\eta_l} s = w^\top B\_{kl}w, \quad
\partial\_{\eta_k\eta_l} w = \Sigma^{-1}B\_{kl}w, \quad
\partial\_{\eta_k\eta_l}\Sigma^{-1} = \Sigma^{-1}B\_{kl}\Sigma^{-1}.\$\$
The same middle matrix serves all three, which is what keeps the
assembly short. \\c\\ and \\d\\ then follow from the quotient rule with
\\\nu+p\\ linear in \\\nu\\ and constant in everything else.
