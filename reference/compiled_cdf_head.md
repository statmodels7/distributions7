# Power-Law Head of a Compiled CDF Derivative

Computes the part of \\\partial_I F(q) = \int_0^q f B_I\\ below the cut
\\y_0\\ that `cdf_rule_cpp()` reports for a family on \\(0, \infty)\\
whose density behaves like \\y^{\alpha - 1}\\ near zero (a gamma or a
generalized gamma with a small shape), where part of the mass lies below
the smallest representable double.

## Usage

``` r
compiled_cdf_head(distrib, rule, th, order, params)
```

## Arguments

- distrib:

  A family named in
  [`compiled_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf_deriv_k.md).

- rule:

  The list that `cdf_rule_cpp()` returns.

- th:

  The named list of parameter vectors, one value per point.

- order:

  The derivative order, 3 or 4.

- params:

  The parameter names.

## Value

`NULL` when no point has a head; otherwise a list with `idx`, the points
that have one, and `integral`, a function of a multi-index that returns
the head's integral at those points.

## Details

On \\(0, y_0\]\\, with \\L = \log y\\, every derivative of \\\ell\\ in
the parameters is linear in \\L\\, so it is read at \\10^{-290}\\ and 15
units of \\L\\ below it and extended to any \\L\\; \\B_I\\ is then a
polynomial in \\L\\ of degree at most four, and
\$\$\int\_{-\infty}^{L_0} e^{c + \alpha u} B_I(u)\\ du = \frac{f(y_0)\\
y_0}{\alpha} \int_0^\infty e^{-t} B_I(L_0 - t/\alpha)\\ dt\$\$ is taken
by the three-point Gauss-Laguerre rule, exact for a polynomial of degree
five. The cut, \\\alpha\\ and \\f(y_0)\\y_0\\ come from
`cdf_rule_cpp()`, which uses the same head at orders one and two.
