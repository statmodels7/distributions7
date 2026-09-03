# The Derivative of the Expected Information

\\\partial\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\, one component
per pair \\(a,b)\\ and differentiating parameter \\c\\.

## Usage

``` r
distrib_dexpected_hessian(
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

  `"parameter"` for \\\partial/\partial\theta_c\\, `"link"` for
  \\\partial/\partial\eta_c\\ of the link-scale expected information.

- approx, nsim:

  Passed to
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md).

- ...:

  Passed to methods.

## Value

A named list of numeric vectors, keyed as
[`dexpected_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md).

## Details

The components are symmetric in \\(a,b)\\ and NOT in \\c\\: writing
\$\$\frac{\partial}{\partial\theta_c}\mathbb{E}\[\ell\_{ab}\] =
\mathbb{E}\[\ell\_{abc}\] + \mathbb{E}\[\ell\_{ab}\ell\_{c}\],\$\$ the
first term is fully symmetric and the second is not, so the result is
keyed by
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
rather than by the sorted triples
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
uses at order three.

**The default method differences the family's own expected
information**, one central stencil per parameter, which is a single
difference of an analytic quantity wherever that quantity is a
written-out formula. It is refused where it is not: see
[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md).

On `scale = "link"` the difference is taken along the free scale of the
parameter being differentiated, and the expected information is read on
the link scale at each of the two points, so the chain rule is never
written out here and cannot disagree with the one
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
already applies.

## See also

[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md)

## Examples

``` r
d <- gaussian1_distrib()
str(distrib_dexpected_hessian(d, 0, list(mu = 0, sigma = 1)))
#> List of 6
#>  $ mu_mu_mu         : num 0
#>  $ mu_mu_sigma      : num 2
#>  $ sigma_sigma_mu   : num 0
#>  $ sigma_sigma_sigma: num 4
#>  $ mu_sigma_mu      : num 0
#>  $ mu_sigma_sigma   : num 0
```
