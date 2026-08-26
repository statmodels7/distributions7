# CDF Derivatives of a Discrete Distribution

Evaluates \\\partial^I F(q)\\ for a discrete family as the finite sum
\$\$\partial^I F(q) = \sum\_{y \le q} f(y)\\ \frac{\partial^I
f}{f}(y).\$\$ Nothing is differenced and nothing is integrated: the sum
is exact.

## Usage

``` r
discrete_cdf_deriv(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib`.

- q:

  A numeric vector of quantiles. Each is handled separately, the support
  being walked from the lower bound up to it, so the cost grows with the
  largest quantile asked for.

- theta:

  A named list of parameters on the parameter scale.

- order:

  The derivative order, 1 or 2. At order 1 the components are the
  parameters and at order 2 they are
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale.
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
puts them on the requested tail.

## Where the identity comes from

The region of integration does not depend on \\\theta\\, so derivative
and integral exchange and \$\$\frac{\partial^I F(q)}{F(q)} =
\mathbb{E}\\\left\[\frac{\partial^I f}{f} \\\middle\|\\ Y \le
q\right\],\$\$ a partial expectation of exactly the quantity the
Bartlett lemma expands: the score at first order, and \\\ell^{(ij)} +
\ell^{(i)}\ell^{(j)}\\ at second. For a discrete family that conditional
expectation is a finite sum whenever the support has a finite lower
bound, which the discrete class requires, so the identity is used as it
stands.

## A trap for whoever writes a test

Checking this against the partial-expectation sum proves nothing: it is
the same sum computed twice. A discrete implementation has to be checked
against finite differences of the cdf, and a continuous one against the
partial expectation, so that the two sides of the comparison share no
arithmetic.

## Notation

\\f\\ is the mass function, \\F\\ the distribution function, \\\ell =
\log f\\ and \\\partial^I\\ a derivative with respect to a multi-index
of parameters.

## See also

[`distrib_grad_cdf.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.discrete_distrib.md),
its caller;
[`numerical_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv.md),
the continuous route;
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).
