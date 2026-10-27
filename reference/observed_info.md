# The Observed Information

Selects the observed information, the negative Hessian of the
log-likelihood at the data, \$\$-\frac{\partial^2
\ell}{\partial\eta\\\partial\eta^\top},\$\$ evaluated at the estimates.
The matrix is computed from the second derivatives of the family, which
are analytic where the family writes them and numerical otherwise (see
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)).
The function takes no argument.

It is the default of the `information` argument of
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
and
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
for a family whose log-likelihood is twice differentiable in every
parameter (see
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)).
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
compares it with
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md).

## Usage

``` r
observed_info()
```

## Value

An S7 object of class `ObservedInfo`, which inherits from
[Information](https://statmodels7.github.io/distributions7/reference/Information.md),
with no properties.

## See also

[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
for the other object;
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
for the comparison;
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)
for the function that takes the object.

## Examples

``` r
observed_info()
#> Observed information

set.seed(1)
d <- gamma1_distrib()
fit <- fit_distrib(d, distrib_rng(d, 300, list(mu = 3, phi = 0.6)))
vcov(fit, information = observed_info())
#>                mu           phi
#> mu   1.710378e-02 -9.363797e-19
#> phi -9.363797e-19  2.065086e-03
```
