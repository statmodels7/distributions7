# The Expected (Fisher) Information

Selects the expected information, the Fisher information
\$\$\mathbb{E}\left\[-\frac{\partial^2
\ell}{\partial\eta\\\partial\eta^\top} \right\],\$\$ evaluated at the
estimates. With `method = NULL`, the default, the closed form is used
where the family writes the expected information out (see
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)),
and
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
otherwise. An object in `method` is used for every family, also for one
that has the closed form.

## Usage

``` r
expected_info(method = NULL)
```

## Arguments

- method:

  `NULL` (the default) or an object that selects how the expectation is
  obtained:

  - `NULL`: the closed form of the family where it has one, and
    [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
    otherwise;

  - [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md):
    the outer product of the observed scores, which costs one gradient;

  - [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md):
    quadrature (or summation over the support) of the observed
    information against the density with `identity = "hessian"`, or of
    minus the outer product of the scores with `identity = "score"`;

  - [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md):
    Monte Carlo, with its own sample size `nsim` and `identity`.

  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
  compares the methods. Any other value signals an error.

## Value

An S7 object of class `ExpectedInfo`, which inherits from
[Information](https://statmodels7.github.io/distributions7/reference/Information.md),
with the property `method`.

## See also

[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
for the other object;
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
for the comparison;
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the methods;
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)
for the function that takes the object.

## Examples

``` r
expected_info()
#> Expected information (closed form where the family has one, otherwise the outer product of the scores)
expected_info(method = expect_mc(nsim = 2000))
#> Expected information (expectation by Monte Carlo (nsim = 2000))

set.seed(1)
d <- gamma1_distrib()
fit <- fit_distrib(d, distrib_rng(d, 300, list(mu = 3, phi = 0.6)))
vcov(fit, information = expected_info())
#>             mu         phi
#> mu  0.01710378 0.000000000
#> phi 0.00000000 0.002065084
```
