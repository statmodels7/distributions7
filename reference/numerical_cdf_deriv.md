# Finite-Difference CDF Derivatives

Differentiates the distribution function with respect to the parameters
by central differences of
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md).
This is the route that a continuous family takes when it registers no
closed form, and it is deliberately chosen there over the exact
partial-expectation identity: for a continuous family that identity is
an integral over a semi-infinite region, while the cdf itself is
analytic for every family in the catalog, so differencing it is both
cheaper and more accurate.

## Usage

``` r
numerical_cdf_deriv(
  distrib,
  q,
  theta,
  order = 1L,
  h_rel = .Machine$double.eps^(1/(order + 2)),
  which = NULL
)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors; the step is then chosen elementwise.

- order:

  The derivative order, 1 (the default) or 2.

- h_rel:

  The relative step. It is a single number that defaults to
  \\\varepsilon^{1/(\mathrm{order}+2)}\\, which balances the stencil's
  truncation against its rounding. A step much smaller than the default
  is worse, because the rounding grows as \\h^{-\mathrm{order}}\\.

- which:

  A character vector naming the components to differentiate, or `NULL`
  (the default) for all of them: parameter names at order 1 and
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
  components at order 2. A family with a closed form for part of the
  order passes the rest here, so that only those cost cdf evaluations.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself and not of
its logarithm. At order 1 the components follow `which`, or the
parameter order when `which` is `NULL`. At order 2 they are always in
the order of
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
restricted to the names in `which`.

## The stencils

At order 1 each component is one central difference, \\\\F(\theta + h) -
F(\theta - h)\\/(2h)\\. At order 2 a diagonal component is the
three-point second difference and an off-diagonal one is the four-point
mixed stencil \\\\F(+,+) - F(+,-) - F(-,+) + F(-,-)\\/(4h_ih_j)\\, which
differences two different variables and is therefore one stencil; the
package never composes two first differences in the same variable.

## The step

The relative step is \\\varepsilon^{1/(k+2)}\\, which is
\\6.1\times10^{-6}\\ at order 1 and \\1.2\times10^{-4}\\ at order 2. It
is multiplied, for each parameter \\j\\, by `pmax(1, abs(theta[[j]]))`.
Because `theta` may be vectorized, this product is a vector with one
step per observation, so that every element of a vectorized parameter
receives a step of the right size.

## Where it is used

The shipped univariate families do not call this function for their
gradient or Hessian: the families whose derivative in a shape parameter
has no elementary form (the gamma, chi-squared, generalized gamma, beta,
von Mises, Student t, pseudo-Huber and skew t families) register the
compiled quadrature of
[`compiled_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md).
The function serves continuous families defined by the user, the
wrappers of continuous families, and the first two tables of
[`cdf_tables()`](https://statmodels7.github.io/distributions7/reference/cdf_tables.md)
where a family has no method of its own.

## Notation

\\F\\ is the distribution function, \\\theta\\ the parameter on its own
scale, \\h\\ the step and \\\varepsilon\\ the machine epsilon.

## See also

[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md),
the exact route for a discrete family;
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md),
which puts the result on the requested tail;
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
and
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md).

## Examples

``` r
d <- gaussian1_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

# The function takes one central difference per parameter.
numerical_cdf_deriv(d, q, th, order = 1)
#> $mu
#> [1] -0.1848768 -0.3278664 -0.1218783
#> 
#> $sigma
#> [1]  0.20028320 -0.05464441 -0.17266092
#> 

# The result is compared with the closed form that the Gaussian registers,
# on the same scale.
exact <- distrib_grad_cdf(d, q, th, log = FALSE)
max(abs(numerical_cdf_deriv(d, q, th, 1)$mu / exact$mu - 1))
#> [1] 6.085332e-11

# Only the requested components are differenced.
names(numerical_cdf_deriv(d, q, th, order = 2, which = "mu_mu"))
#> [1] "mu_mu"
```
