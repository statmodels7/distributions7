# Integration Against the Density as the Expected Derivative

Selects numerical integration over the response as the method that
obtains an expected derivative. The quantity integrated against the
density \\f\\ is the one that `identity` names:
\$\$\mathbb{E}\[\partial^k \ell\] = \int \partial^k \ell(y;\theta)\\
f(y;\theta)\\dy \quad \text{for } \texttt{identity = "hessian"},\$\$ the
observed derivative of the requested order, or \$\$-\sum\_{\pi} \int
\prod\_{B \in \pi} \ell_B(y;\theta)\\f(y;\theta)\\dy \quad \text{for }
\texttt{identity = "score"},\$\$ where \\\pi\\ runs over the partitions
of the index set into at least two blocks and \\\ell_B\\ is the
derivative of the log-density with respect to the parameters indexed by
the block \\B\\ (the Bartlett partition identity; at order 2 it reads
\\-\int \ell_i \ell_j f\\dy\\). The integrals are computed by batched
adaptive quadrature for a continuous family and by series summation for
a discrete one. The computation is deterministic. The accuracy of the
quadrature and of the series is set by the arguments of the constructor,
whose defaults are those of
[`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html)
and
[`numericals7::series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.html).

With `identity = "hessian"` the integral needs the derivative of the
target order and costs one integral per component; it is normally the
most accurate method when that derivative is available in closed form.
With `identity = "score"` the target order is never needed, only the
lower ones, and the cost is one integral per partition and component.
The score identity is the one that remains valid when the log-likelihood
has a kink in a parameter: at order 2, for a family with a non-smooth
parameter (see
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md)),
the observed second derivative lacks the point mass at the kink, and
`identity = "hessian"` signals an error that names `identity = "score"`.

With the defaults this is the method that
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
use for an expected derivative when `method` is `NULL`.

## Usage

``` r
expect_integrate(
  identity = c("hessian", "score"),
  quad_atol = 1e-10,
  quad_rtol = 1e-08,
  max_panels = 4096L,
  series_atol = 1e-12,
  series_rtol = 1e-10,
  max_terms = 100000L
)
```

## Arguments

- identity:

  One of `"hessian"` (the default) and `"score"`. See the description
  and
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- quad_atol, quad_rtol:

  The absolute and the relative tolerance of the quadrature of a
  continuous family, single positive finite numbers. A panel of the
  integration is accepted when its error estimate is below
  `quad_atol + quad_rtol * abs(integral)`. The defaults are 1e-10 and
  1e-8.

- max_panels:

  The largest number of panels that the quadrature of one parameter
  configuration may reach, a single positive finite number. A row that
  needs more panels is not accepted by the batched quadrature and is
  integrated by the scalar integrator of stats. The default is 4096.

- series_atol, series_rtol:

  The absolute and the relative tolerance of the series summation of a
  discrete family, single positive finite numbers. The defaults are
  1e-12 and 1e-10.

- max_terms:

  The largest number of terms of a series, a single positive finite
  number. A series that has not met its tolerance after that many terms
  signals an error. The default is 100000.

## Value

An S7 object of class `ExpectIntegrate`, which inherits from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
with the properties `identity`, `quad_atol`, `quad_rtol`, `max_panels`,
`series_atol`, `series_rtol` and `max_terms`.

## See also

[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
for the other methods;
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the comparison;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the functions that take the object.

## Examples

``` r
Dens <- S7::new_class("DensIntegrate", parent = continuous_distrib)
S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
  stats::dnorm(y, theta[[1]], theta[[2]], log = log)
}
d <- Dens(distrib_name = "dens", dimension = "univariate",
          bounds = c(-Inf, Inf), params = c("mu", "sigma"),
          params_interpretation = c(mu = "location", sigma = "scale"),
          n_params = 2L,
          params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
          link_params = list(mu = linkfunctions7::identity_link(),
                             sigma = linkfunctions7::log_link()),
          params_smooth = c(mu = TRUE, sigma = TRUE))
th <- list(mu = 0, sigma = 2)
# The expected information of sigma is 2 / sigma^2 = 0.5 in absolute value,
# by either identity.
distrib_expected_hessian(d, 0, th, method = expect_integrate())$sigma_sigma
#> [1] -0.5
distrib_expected_hessian(d, 0, th,
                         method = expect_integrate(identity = "score"))$sigma_sigma
#> [1] -0.5

# Tighter tolerances for the quadrature are set in the object.
a <- expect_integrate(quad_atol = 1e-12, quad_rtol = 1e-10)
a
#> Expectation by integration against the density (quad_atol = 0.000000000001, quad_rtol = 0.0000000001)
distrib_expected_hessian(d, 0, th, method = a)$sigma_sigma
#> [1] -0.5
```
