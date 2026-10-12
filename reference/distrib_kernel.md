# A Resolved Kernel for One Parameter's Link-Scale Derivatives

The log-density, the score, the curvature, the expected information and
its derivative for one parameter on its unconstrained scale, as five
functions with everything that does not depend on the data already
resolved.

## Usage

``` r
distrib_kernel(distrib, param)
```

## Arguments

- distrib:

  A univariate distribution object.

- param:

  The name of the parameter with respect to whose unconstrained scale
  the derivatives are taken.

## Value

A list of five functions of `(y, theta, eta)`: `logdens`, `score`,
`curvature`, `information` and `dinformation`.

## Details

The generic route validates its arguments, aligns `theta` by name,
dispatches, and assembles every component of the requested order. A
score-driven filter evaluates the score once per observation, at a
predictor that it has just produced, so these steps cannot be vectorized
away and their cost dominates.

This resolves the family's methods and the link's once, and applies the
chain rule for the single component wanted rather than for all of them.
The Jacobian of the parametrization is diagonal, so with \\\theta_p =
h(\eta_p)\\,

\$\$\frac{\partial \ell}{\partial \eta_p} = \ell_p\\h'(\eta_p), \qquad
\frac{\partial^2 \ell}{\partial \eta_p^2} = \ell\_{pp}\\h'(\eta_p)^2 +
\ell_p\\h''(\eta_p),\$\$

which is what
[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
computes for those two components and nothing else.

The caller takes on the checks that the generic would have made. `theta`
must already be a list named by the parameters of the family, in the
family's own order, with values of a length that the family accepts
against `y`; nothing is checked. The entry named `param` is replaced, so
its value on the way in is immaterial. A list without names is not
replaced at that entry, and a list whose names are out of order is read
by position. The inverse link is clamped strictly inside its bounds
exactly as
[`linkfunctions7::linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.html)
does, which correctness requires.

The kernel also carries the expected information of the same coordinate
and its derivative in that coordinate, which a filter whose score is
scaled by a power of the information reads at each step. On the link
scale the information is \\I_p =
-\mathbb{E}\[\ell\_{pp}\]\\h'(\eta_p)^2\\, with no term in \\h''\\
because the score has mean zero, and

\$\$\frac{\partial I_p}{\partial \eta_p} =
-\frac{\partial\mathbb{E}\[\ell\_{pp}\]}{\partial\theta_p}\\h'^3 -
2\\\mathbb{E}\[\ell\_{pp}\]\\h'\\h''.\$\$

The two methods are looked up when first called, so a family without an
expected information still yields the other three functions.

## See also

[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md),
[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)

## Examples

``` r
d <- gaussian1_distrib()
k <- distrib_kernel(d, "sigma")
th <- list(mu = 0, sigma = 1)
k$score(0.7, th, log(1.4))
#> [1] -0.75
# the same number the generic gives
distrib_gradient(d, 0.7, list(mu = 0, sigma = 1.4),
                 scale = "link")[["sigma"]]
#> [1] -0.75
```
