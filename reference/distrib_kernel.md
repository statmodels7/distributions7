# A Resolved Kernel for One Parameter's Link-Scale Derivatives

The log-density, the score and the curvature in ONE parameter's
unconstrained scale, as three functions with everything that does not
depend on the data already resolved.

## Usage

``` r
distrib_kernel(distrib, param)
```

## Arguments

- distrib:

  A univariate distribution object.

- param:

  The name of the parameter whose unconstrained scale the derivatives
  are taken with respect to.

## Value

A list of three functions of `(y, theta, eta)`: `logdens`, `score` and
`curvature`.

## Details

The generic route is the right one almost everywhere: it validates its
arguments, aligns `theta` by name, dispatches, and assembles every
component of the requested order. A recursion that calls back once per
observation cannot afford any of that. A score-driven filter evaluates
the score at a predictor it has just produced, so the call cannot be
vectorized away, and profiling a fitted model put S7 dispatch, the
argument checking and the name arithmetic at the whole of its cost.

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

The bargain is that the caller takes on what the generic was doing.
`theta` must already be a list in the family's own order, its values
unnamed and of a length the family accepts against `y`; nothing is
checked. The entry for `param` is replaced, so its value on the way in
is immaterial. The inverse link is clamped strictly inside its bounds
exactly as
[`linkfunctions7::linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.html)
does, because that is a correctness property and not an optimization.

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
