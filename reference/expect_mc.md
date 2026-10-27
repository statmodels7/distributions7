# Monte Carlo as the Expected Derivative

Selects Monte Carlo as the method that obtains an expected derivative.
The function draws `nsim` observations from the distribution at each
parameter configuration, with
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
and averages over them the quantity that `identity` names:
\$\$\widehat{\mathbb{E}}\[\partial^k \ell\] = \frac{1}{n\_{sim}}
\sum\_{m=1}^{n\_{sim}} \partial^k \ell(y_m;\theta)\$\$ for
`identity = "hessian"`, the observed derivative of the requested order,
and the same average of \\-\sum\_{\pi} \prod\_{B \in \pi}
\ell_B(y_m;\theta)\\ for `identity = "score"`, the Bartlett partition
identity (see
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)).
Both estimate the same quantity as
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
with the same `identity`. The error decreases as \\1/\sqrt{n\_{sim}}\\,
and the result is reproducible only if the seed is fixed. The cost is
that of
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
times `nsim`, per distinct parameter configuration. The choice is
adequate where quadrature performs poorly, as with heavy tails or an
awkward support.

For a family with a non-smooth parameter (see
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md))
the average of the observed second derivative omits the point mass at
the kink, so that `identity = "hessian"` at order 2 signals an error
that names `identity = "score"`.

## Usage

``` r
expect_mc(nsim = 10000, identity = c("hessian", "score"))
```

## Arguments

- nsim:

  The number of draws per parameter configuration, a single positive
  finite number. The default is 10000. A vector, a missing or non-finite
  value, or a number below 1 signals an error.

- identity:

  One of `"hessian"` (the default) and `"score"`. See the description
  and
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

## Value

An S7 object of class `ExpectMC`, which inherits from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
with the properties `nsim` and `identity`.

## See also

[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
and
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
for the other methods;
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the comparison;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
and
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for the functions that take the object.

## Examples

``` r
a <- expect_mc(nsim = 2000)
a
#> Expectation by Monte Carlo (nsim = 2000)
expect_mc(nsim = 2000, identity = "score")
#> Expectation by Monte Carlo (nsim = 2000, identity = score)

# A sample is drawn at every call, so the seed fixes the result.
d <- gaussian1_distrib()
th <- list(mu = 0, sigma = 2)
set.seed(1)
distrib_deriv3(d, 0, th, expected = TRUE, method = a)$mu_mu_sigma
#> [1] 0.25

# The step of Fisher scoring can use the same object.
fisher_scoring(method = expect_mc(nsim = 1000))
#> Fisher scoring
#>   expected information: by Monte Carlo (nsim = 1000)
```
