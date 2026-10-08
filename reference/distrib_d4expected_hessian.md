# The Fourth Derivative of the Expected Information

\\\partial^4\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\\partial\theta_d\\
\partial\theta_e\\\partial\theta_f\\, one component per pair \\(a,b)\\
and per unordered quadruple \\(c,d,e,f)\\.

## Usage

``` r
distrib_d4expected_hessian(
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
[`d4expected_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

## Details

The components are symmetric in \\(a,b)\\ and in \\(c,d,e,f)\\
separately, and are keyed by
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md).

The default method applies one second-order central stencil to
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md):
the three-point stencil along one parameter for a repeated pair
\\(e,e)\\, the four-point mixed stencil for two distinct parameters
\\(e,f)\\. It requires a family that registers an analytic second
derivative of its expected information and signals an error otherwise.
The step is the fourth root of machine epsilon relative to the
coordinate, kept inside a finite bound on the parameter scale by
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md).

## See also

[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md),
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md),
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)

## Examples

``` r
d <- poisson_distrib()
str(distrib_d4expected_hessian(d, 2, list(mu = 3), scale = "link"))
#> List of 1
#>  $ mu_mu_mu_mu_mu_mu: num -3
```
