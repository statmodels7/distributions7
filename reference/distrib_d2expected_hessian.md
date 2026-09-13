# The Second Derivative of the Expected Information

\\\partial^2\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\\partial\theta_d\\,
one component per pair \\(a,b)\\ and per pair \\(c,d)\\.

## Usage

``` r
distrib_d2expected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  approx = c("opg", "bartlett", "integrate", "mc"),
  nsim = 10000,
  ...
)
```

## Arguments

- distrib:

  A distribution object inheriting from `distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, each of length 1 or `length(y)`.

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Accepted for symmetry with
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md);
  no method reads them.

- ...:

  Passed to methods.

## Value

A named list of numeric vectors, keyed as
[`d2expected_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

Differentiating the identity of
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
once more moves the measure a second time,
\$\$\frac{\partial^2}{\partial\theta_c\\\partial\theta_d}\mathbb{E}\[\ell\_{ab}\]
= \mathbb{E}\[\ell\_{abcd}\] + \mathbb{E}\[\ell\_{abd}\ell\_{c}\] +
\mathbb{E}\[\ell\_{abc}\ell\_{d}\] +
\mathbb{E}\[\ell\_{ab}\ell\_{cd}\] +
\mathbb{E}\[\ell\_{ab}\ell\_{c}\ell\_{d}\],\$\$ and no Bartlett identity
isolates those moments. A family supplies the components as ordinary
derivatives of its own written-out expected information; the components
are symmetric in \\(a,b)\\ and in \\(c,d)\\ separately, and are keyed by
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

**There is no numerical default.** A difference of
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
whose own default is already a difference, would be the nested
differencing the package forbids, so the base method signals an error
and a family that does not register one has no second derivative. The
families that do are
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md),
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md),
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md),
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
and
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md),
each from a compiled kernel.

On `scale = "link"` the expected information is \\F\_{ab} =
\mathbb{E}\[\ell\_{ab}\]\\h_a' h_b'\\, with no term in \\h''\\ because
\\\mathbb{E}\[\ell_a\] = 0\\, and its derivatives follow by Leibniz's
rule, written once in
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md).

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)

## Examples

``` r
d <- gaussian1_distrib()
str(distrib_d2expected_hessian(d, 0, list(mu = 0, sigma = 1)))
#> List of 9
#>  $ mu_mu_mu_mu            : num 0
#>  $ mu_mu_sigma_sigma      : num -6
#>  $ mu_mu_mu_sigma         : num 0
#>  $ sigma_sigma_mu_mu      : num 0
#>  $ sigma_sigma_sigma_sigma: num -12
#>  $ sigma_sigma_mu_sigma   : num 0
#>  $ mu_sigma_mu_mu         : num 0
#>  $ mu_sigma_sigma_sigma   : num 0
#>  $ mu_sigma_mu_sigma      : num 0
```
