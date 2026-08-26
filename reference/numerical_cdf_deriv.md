# Finite-Difference CDF Derivatives

Differentiates the distribution function with respect to the parameters
by central differences of
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md).
This is the route a continuous family takes when it registers no closed
form, and it is deliberately preferred there to the exact
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

  The relative step. A single number, defaulting to
  \\\varepsilon^{1/(\mathrm{order}+2)}\\, which balances the stencil's
  truncation against its rounding. A step much smaller than the default
  is worse, the rounding growing as \\h^{-\mathrm{order}}\\.

- which:

  A character vector naming the components to differentiate, or `NULL`
  (the default) for all of them: parameter names at order 1 and
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
  components at order 2. A family with a closed form for part of the
  order passes the rest here, so that only those cost cdf evaluations.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself and not of
its logarithm, in the order `which` gives or in the enumeration order
when it is `NULL`.

## The stencils

At order 1 each component is one central difference, \\\\F(\theta + h) -
F(\theta - h)\\/(2h)\\. At order 2 a diagonal component is the
three-point second difference and an off-diagonal one is the four-point
mixed stencil \\\\F(+,+) - F(+,-) - F(-,+) + F(-,-)\\/(4h_ih_j)\\, which
differences two different variables and is therefore one stencil; the
package never composes two first differences in the same variable.

## The step

The relative step is \\\varepsilon^{1/(k+2)}\\, which is
\\6.1\times10^{-6}\\ at order 1 and \\1.2\times10^{-4}\\ at order 2, and
it is scaled per component by `pmax(1, abs(theta[[j]]))`. **One step is
chosen per observation and not per parameter**: `theta` may be
vectorized, and a step read off its first element would be the wrong
size everywhere else.

## What it costs and what it delivers

Measured on a Gaussian at 1000 quantiles, against the closed form the
family registers: the relative error is \\6.1\times10^{-11}\\ at order 1
and \\1.7\times10^{-7}\\ at order 2, and the calls take 0.29 ms and 0.71
ms against the closed form's 0.16 ms. Of the 42 univariate families, 8
reach this function for their gradient (beta1, beta2, chisq, gamma1,
gamma2, gengamma1 and the two von Mises), where the derivative of an
incomplete gamma or beta in its shape is hypergeometric or the cdf is
itself a quadrature.

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

# One central difference per parameter.
numerical_cdf_deriv(d, q, th, order = 1)
#> $mu
#> [1] -0.1848768 -0.3278664 -0.1218783
#> 
#> $sigma
#> [1]  0.20028320 -0.05464441 -0.17266092
#> 

# Against the closed form the Gaussian registers, on the same scale.
exact <- distrib_grad_cdf(d, q, th, log = FALSE)
max(abs(numerical_cdf_deriv(d, q, th, 1)$mu / exact$mu - 1))
#> [1] 6.085332e-11

# Only the components asked for are differenced.
names(numerical_cdf_deriv(d, q, th, order = 2, which = "mu_mu"))
#> [1] "mu_mu"
```
