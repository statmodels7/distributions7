# The Weights That Carry the Tail of the Skew t Past the Rule

Computes, for each degrees of freedom \\\nu\\, the weights of the nine
anchor nodes of
[`loc_scale_tail_anchors()`](https://statmodels7.github.io/distributions7/reference/loc_scale_tail_anchors.md)
in each tail, so that the sum over the nodes includes the nodes the rule
does not have.

## Usage

``` r
loc_scale_tail_weights(nu)
```

## Arguments

- nu:

  A numeric vector of degrees of freedom.

## Value

A matrix with nine rows and one column for each element of `nu`.

## Details

Past \\\|z\| = 10^{60}\\ the density of the skew t is \\f(z) =
f(z_i)(\|z\|/z_i)^{-(\nu+1)}\\ to a relative \\10^{-26}\\ or better from
\\z_i = 10^{13}\\ on, and each derivative of the log-density is constant
except the derivative in \\\nu\\, which is linear in \\u = \log\|z\|\\.
Every integrand of the expected information and of its derivatives is
therefore \\f\\ times a polynomial \\q(u)\\ of degree at most six. The
sum of \\w_k f(z_k) q(u_k)\\ over the nodes \\k \> k\_{\max}\\ that the
rule omits is then, with \\q\\ written in the Lagrange basis \\L_i\\ on
the anchors \\u_i = \log z_i\\, \$\$\sum_i f(z_i)\\ q(u_i)\\
\Omega_i,\qquad \Omega_i = z_i \sum\_{k \> k\_{\max}} h c_0 \cosh(t_k)\\
e^{-\nu (u_k - u_i)} L_i(u_k),\$\$ with \\u_k = c_0 \sinh t_k\\ the grid
of the rule, so that \\\Omega_i\\ depends on \\\nu\\ alone. The sum over
\\k\\ stops where \\\nu (u_k - u\_{k\_{\max}})\\ reaches 100. At \\\nu
\ge 1\\ the tail is below \\10^{-45}\\ of the integral and the weights
are zero.
