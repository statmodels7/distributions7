# The Expected Information of a Folded Family

Computes the expected information of a
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
family (`order = 0`), its first derivatives (`order = 1`) or its first
and second derivatives (`order = 2`), on the parameter scale, as
expectations under the folded law.

## Usage

``` r
folded_expected(distrib, y, theta, order)

fold_kink_parts(distrib, th, center, order)
```

## Arguments

- distrib:

  A
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
  family.

- y:

  The response, read for its length.

- theta:

  The parameters.

- order:

  `0L` for the expected information, `1L` for its first derivatives,
  `2L` for its first and second derivatives.

- th:

  One row of parameters, in the order of the family's.

- center:

  The parent's center at that row.

## Value

A named list keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
(`order = 0`),
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
(`order = 1`), or
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
followed by
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
(`order = 2`), each component of length `length(y)`.

## Details

An expectation under the folded law is an integral over \\y \> 0\\
against the folded density \\f(y) + f(-y)\\. The integrand turns sharply
near zero, where the weight \\f(y)/(f(y)+f(-y))\\ goes from one half to
one over a width of order \\s^2/\|c\|\\, \\c\\ and \\s\\ a center and a
scale of the parent read from the scalar registry. The integral is taken
over \\y \> 0\\ by a tanh-sinh rule on \\\[0, \|c\|\]\\, whose nodes
cluster at both ends, and the rule of
[`loc_scale_rule()`](https://statmodels7.github.io/distributions7/reference/loc_scale_rule.md)
scaled by \\s\\ on \\\[\|c\|, \infty)\\; the split at \\\|c\|\\ also
takes in a parent with a kink at its center. The nodes are the scalar
registry's. The expectations are formed through the second Bartlett
identity, \$\$E\[\ell\_{ij}\] = -E\[\ell_i \ell_j\], \qquad \partial_c
E\[\ell\_{ij}\] = -E\[\ell\_{ic}\ell_j + \ell_i\ell\_{jc} +
\ell_i\ell_j\ell_c\],\$\$ so that the folded score and Hessian are all
it reads at orders 0 and 1. The second derivatives differentiate
\\-E\[G\]\\, \\G = \ell_i\ell_j\\, twice under the integral,
\$\$\partial\_{cd} E\[G\] = E\[G\_{cd} + G_c\ell_d + G_d\ell_c +
G(\ell_c\ell_d + \ell\_{cd})\],\$\$ and read the folded third
derivatives as well. At order 2 a parent without a route signals an
error. The sums are accumulated as the scalar registry accumulates them.
A parent that the registry does not cover, or that is not a
location-scale family on the real line, is integrated by
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
with `approx = "integrate"` instead.

A parent with a kink at its center \\c\\ (a parameter that is not
smooth, as in
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md),
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
and
[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md))
gives a folded score that jumps at \\y^\* = \|c\|\\, a point that moves
with \\c\\. With \\F = L\\G\\ the integrand, \\J(X) = X(y^{\*-}) -
X(y^{\*+})\\ and \\a_k = \partial y^\*/\partial \theta_k =
\mathrm{sign}(c)\\ for the center and zero otherwise, \$\$\partial_k
\int F = \int \partial_k F + J(F)\\a_k,\$\$ \$\$\partial\_{kl} \int F =
\int \partial\_{kl} F + J(\partial_k F)\\a_l + J(\partial_l F)\\a_k +
J(\partial_y F)\\a_k a_l,\$\$ the one-sided values read four units in
the last place either side of \\y^\*\\. `fold_kink_parts()` returns
those values: the folded density at \\y^\*\\, the folded score and
Hessian, and, at order 2, the response derivative \\\ell_y\\ and the
mixed derivatives \\\ell\_{iy}\\ from the parent's at the two preimages,
or `NULL` where the parent has no kink or \\c = 0\\.

## See also

[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
