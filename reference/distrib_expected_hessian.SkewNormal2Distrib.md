# Skew Normal Expected Information in the Centered Parametrization

The expected information and its first two derivatives in the
parameters. For \\\|\gamma_1\| \<\\
[`sn2_ge()`](https://statmodels7.github.io/distributions7/reference/sn2_ge.md)
they come from their series in \\r = (\gamma_1/c)^{1/3}\\; elsewhere
from the quadrature of
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
over the family's own observed derivatives, the family being
location-scale in \\(\mu, \sigma)\\ at fixed \\\gamma_1\\.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object.

- y:

  A numeric vector of observations, read for its length.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`.

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Accepted for the generics' signatures and unused.

- ...:

  Unused.

- threads:

  A single positive integer, passed to the kernels behind the
  quadrature.

## Value

A named list keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

The series is the series in \\r\\ of the observed components integrated
term by term against the series of the density, with gaussian moments,
computed offline at 60 digits. It is asymptotic rather than convergent,
and below the bound its terms fall under \\10^{-20}\\ before they turn.
The information is finite at \\\gamma_1 = 0\\, where it is
\\\mathrm{diag}(1, 2, 1/6)/\sigma^2\\ with the sign of a Hessian, but it
is not analytic in \\\gamma_1\\ there: \\E\[\ell\_{\gamma_1\gamma_1}\]\\
carries a term in \\\gamma_1^{2/3}\\ and \\E\[\ell\_{\mu\gamma_1}\]\\
one in \\\gamma_1^{4/3}\\. Its derivatives in \\\gamma_1\\ are therefore
infinite at zero skewness, where the two derivative methods signal an
error.

## See also

[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md),
[`sn2_ge()`](https://statmodels7.github.io/distributions7/reference/sn2_ge.md).

## Examples

``` r
d <- skewnormal2_distrib()
# finite at zero skewness, where the direct parametrization's is singular
distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, gamma1 = 0))
#> $mu_mu
#> [1] -1
#> 
#> $sigma_sigma
#> [1] -2
#> 
#> $gamma1_gamma1
#> [1] -0.1666667
#> 
#> $mu_sigma
#> [1] 0
#> 
#> $mu_gamma1
#> [1] 1.166815e-61
#> 
#> $sigma_gamma1
#> [1] 0
#> 
```
