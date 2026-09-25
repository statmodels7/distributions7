# The Expected Information of a Location-Scale Family

Computes the expected information of a family whose first two parameters
are a location \\\mu\\ and a scale \\\sigma\\, and optionally its first
and second derivatives in the parameters, by one quadrature per distinct
value of the remaining shape parameters.

## Usage

``` r
loc_scale_expected(distrib, theta, order, n, threads = 1L)
```

## Arguments

- distrib:

  A family whose first two parameters are a location and a scale, and
  which carries analytic derivatives to the order needed.

- theta:

  The aligned parameter list, each component of the common length `n`.

- order:

  `0L` for the expected information, `1L` for its first derivatives,
  `2L` for its second. At order 2 the first derivatives are returned
  beside the second, from the same evaluation.

- n:

  The number of observations.

- threads:

  The thread count passed to the family's kernels.

## Value

A named list: at order 0 keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
at order 1 as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
at order 2 as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
and
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
together. Each component has length `n` and is on the parameter scale.

## Details

When the log-density is \\\ell(y;\mu,\sigma,s) = -\log\sigma + g(z;s)\\
with \\z = (y-\mu)/\sigma\\, every derivative of \\\ell\\ carries one
factor \\1/\sigma\\ per index on \\\mu\\ or \\\sigma\\ and is otherwise
a function of \\z\\ and \\s\\ alone. The same holds after taking the
expectation, so for any component
\$\$\partial^{I}\\\mathbb{E}\[\ell\_{ab}\](\mu,\sigma,s) =
\sigma^{-k}\\\partial^{I}\\\mathbb{E}\[\ell\_{ab}\](0,1,s),\$\$ with
\\k\\ the number of indices among \\a\\, \\b\\ and \\I\\ that fall on
\\\mu\\ or \\\sigma\\. The quantity at \\(0,1,s)\\ does not depend on
the location or on the scale, so it is computed once for each distinct
shape and not once for each observation.

At \\(0,1,s)\\ each quantity is an integral over \\z\\ of the family's
own analytic derivatives against its density \\f\\. Differentiating
under the integral moves the measure as well as the integrand:
\$\$\mathbb{E}\[\ell\_{ab}\] = \int \ell\_{ab}\\f,\$\$
\$\$\partial_c\\\mathbb{E}\[\ell\_{ab}\] = \int (\ell\_{abc} +
\ell\_{ab}\ell_c)\\f,\$\$ \$\$\partial\_{cd}\\\mathbb{E}\[\ell\_{ab}\] =
\int (\ell\_{abcd} + \ell\_{abc}\ell_d + \ell\_{abd}\ell_c +
\ell\_{ab}\ell\_{cd} + \ell\_{ab}\ell_c\ell_d)\\f.\$\$

The integrals are taken by an exp-sinh rule on each half-line, from
[`loc_scale_rule()`](https://statmodels7.github.io/distributions7/reference/loc_scale_rule.md).
Its nodes cluster double-exponentially at \\z = 0\\, which is where the
three families that use it put their one sharp feature (the step of the
skewing function at a large shape, and the core of the pseudo-Huber at a
small one), and the same transformation turns an algebraic tail into an
exponential one, which is what the Student t's heavy tail needs.

A node where the density is exactly zero contributes nothing, whatever
the derivatives read there. Any other non-finite contribution makes the
whole shape row `NA` rather than a sum over the finite part.

## See also

[`loc_scale_rule()`](https://statmodels7.github.io/distributions7/reference/loc_scale_rule.md)
for the quadrature,
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
for the generics this serves.
