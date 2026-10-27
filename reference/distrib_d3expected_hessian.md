# The Third Derivative of the Expected Information

Returns \\\partial^3\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\
\partial\theta_d\\\partial\theta_e\\, with one component per pair
\\(a,b)\\ and per unordered triple \\(c,d,e)\\.

## Usage

``` r
distrib_d3expected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  method = NULL,
  ...
)
```

## Arguments

- distrib:

  A distribution object inheriting from `distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters; each has length 1 or `length(y)`.

- scale:

  `"parameter"` or `"link"`.

- method:

  `NULL` (the default) or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
  as for
  [`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md).
  An object returns the derivative of the matrix that the method
  defines, by Leibniz's formula for
  [`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
  and with the moved density for
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
  and
  [`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md);
  with `identity = "hessian"` the fourth derivative needs observed
  derivatives of order 6 and signals an error that names
  `identity = "score"`. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Passed to methods.

## Value

A named list of numeric vectors, keyed as
[`d3expected_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md).

## Details

The components are symmetric in \\(a,b)\\ and in \\(c,d,e)\\ separately,
and are keyed by
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md).

The univariate families shipped with the package compute the components
in closed form, in compiled code. The default method, which serves the
wrappers and families defined by the user, takes one central difference
of
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
along each parameter, so it requires a family that registers an analytic
second derivative of its expected information and signals an error
otherwise. The step is the cube root of machine epsilon relative to the
coordinate, kept inside a finite bound on the parameter scale by
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md).

On `scale = "link"` the difference of the link-scale second derivative
is taken along the free coordinate of the parameter, so the chain rule
is the one that
[`dexpected_link()`](https://statmodels7.github.io/distributions7/reference/dexpected_link.md)
already applies.

## See also

[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md),
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md),
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)

## Examples

``` r
d <- gaussian1_distrib()
str(distrib_d3expected_hessian(d, 0, list(mu = 0, sigma = 1)))
#> List of 12
#>  $ mu_mu_mu_mu_mu               : num 0
#>  $ mu_mu_mu_mu_sigma            : num 0
#>  $ mu_mu_mu_sigma_sigma         : num 0
#>  $ mu_mu_sigma_sigma_sigma      : num 24
#>  $ sigma_sigma_mu_mu_mu         : num 0
#>  $ sigma_sigma_mu_mu_sigma      : num 0
#>  $ sigma_sigma_mu_sigma_sigma   : num 0
#>  $ sigma_sigma_sigma_sigma_sigma: num 48
#>  $ mu_sigma_mu_mu_mu            : num 0
#>  $ mu_sigma_mu_mu_sigma         : num 0
#>  $ mu_sigma_mu_sigma_sigma      : num 0
#>  $ mu_sigma_sigma_sigma_sigma   : num 0
```
