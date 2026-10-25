# The Third or Fourth Derivative of a Location-Scale Expected Information

Computes the components of
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md)
or
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md)
on the parameter scale for a family whose first two parameters are a
location and a scale, by the quadrature of
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
applied to \\E\_{ab} = -E\[\ell_a \ell_b\]\\.

## Usage

``` r
loc_scale_d34expected(distrib, theta, order, n, threads = 1L)
```

## Arguments

- distrib:

  A distribution object of a location-scale family.

- theta:

  A named list of parameters.

- order:

  `3L` or `4L`.

- n:

  The number of observations.

- threads:

  The thread count passed to the family's derivatives.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

## Details

Differentiating under the integral, \\\partial_S E\[G\] = E\[\sum\_{T
\subseteq S} \partial_T G \\ B\_{S \setminus T}\]\\ with \\B_U = f^{-1}
\partial_U f\\, a sum over the set partitions of \\U\\. The form
\\-E\[\ell_a \ell_b\]\\ reads the derivatives of the log-density up to
order five at the fourth derivative, where \\E\[\ell\_{ab}\]\\ would
read order six. A component with an index on the location is zero, and
\\j\\ indices on the scale contribute the factor \\(-1)^j (k)\_j
\sigma^{-(k+j)}\\, with \\(k)\_j\\ the rising factorial and \\k\\ the
number of indices among \\a, b\\ on the location or the scale. The
quadrature runs once per distinct value of the shapes.
