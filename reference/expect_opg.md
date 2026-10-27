# The Outer Product of the Scores as the Expected Information

Selects the empirical outer product of the observed scores (OPG, the
estimator of Berndt, Hall, Hall and Hausman) as the estimate of an
expected derivative. No expectation is taken and nothing is simulated:
the quantity that the score identity (see
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md))
averages over the model is read at the observations in hand. At order 2
the expected Hessian is replaced, observation by observation, by
\$\$-\\\ell_i(y)\\\ell_j(y),\$\$ where \\\ell_i\\ is the derivative of
the log-density with respect to the \\i\\th parameter. At order 3 the
expected third derivative is replaced by \$\$-\left(\ell\_{ij}\ell_k +
\ell\_{ik}\ell_j + \ell\_{jk}\ell_i + \ell_i\ell_j\ell_k\right),\$\$ and
at order 4 by the corresponding sum over the partitions of the four
indices into at least two blocks, with the sign changed. Every term is a
product of observed derivatives of lower order than the target, so the
cost is one evaluation of each of those orders and the target order
itself is never needed.

The matrix at order 2 is positive semidefinite and has the same fixed
point as the expected information, which makes it adequate for the step
of Fisher scoring. As an estimate of the information it has the variance
of a single product of scores per observation, and at orders 3 and 4 the
variance of a product of up to four derivatives; the sum over \\n\\
observations converges to \\n\\ times the expectation at the rate
\\n^{-1/2}\\. A standard error should be read from the observed
information (see
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)).

This is the default of
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for a family without a closed form, and of
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md).
At orders 3 and 4 it is never the default and is never replaced by
another method: a call with `method = expect_opg()` returns the
empirical sums of products described above.
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
compares it with the other methods.

## Usage

``` r
expect_opg()
```

## Value

An S7 object of class `ExpectOPG`, which inherits from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
with no properties.

## See also

[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
for the other methods;
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the comparison;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
and
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for the functions that take the object.

## Examples

``` r
# A family defined by its density alone has no closed-form information.
Dens <- S7::new_class("DensOPG", parent = continuous_distrib)
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
distrib_expected_hessian(d, c(-1, 0.5, 2), th, method = expect_opg())$mu_mu
#> [1] -0.062500 -0.015625 -0.250000

# The same object is accepted at the higher orders.
distrib_deriv3(d, c(-1, 0.5, 2), th, expected = TRUE,
               method = expect_opg())$mu_mu_sigma
#> [1]  0.05468750 -0.07861328  0.50000001

# And it selects the step of Fisher scoring.
fisher_scoring(method = expect_opg())
#> Fisher scoring
#>   expected information: by the outer product of the scores
```
