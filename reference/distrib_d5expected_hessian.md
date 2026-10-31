# The Fifth Derivative of the Expected Information

Returns \\\partial^5\\\mathbb{E}\[\ell\_{ab}\]/\partial\theta_c\\
\partial\theta_d\\\partial\theta_e\\\partial\theta_f\\\partial\theta_g\\,
with one component per pair \\(a,b)\\ and per unordered quintuple
\\(c,d,e,f,g)\\.

## Usage

``` r
distrib_d5expected_hessian(
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
  [`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md).
  See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Passed to methods; the default method accepts `accuracy`, the order of
  accuracy of the stencil passed to
  [`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html)
  (4 by default).

## Value

A named list of numeric vectors, keyed as
[`d5expected_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/d5expected_names.md).

## Details

The components are symmetric in \\(a,b)\\ and in \\(c,d,e,f,g)\\
separately, and are keyed by
[`d5expected_names()`](https://statmodels7.github.io/distributions7/reference/d5expected_names.md).

No family computes this order in closed form. The default method applies
one tensor-product central stencil of
[`numericals7::fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.html)
to the highest derivative of the expected information that the family
computes itself: a first difference of
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md)
for the univariate families shipped with the package, which compute the
fourth order in closed form, and for a wrapper of such a family; a
stencil of order two on
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md)
or of order three on
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
for a family that registers only those. A wrapper never differences from
an order above the one of the family that it wraps, since its own orders
are built from the parent's. A family without an analytic second
derivative of its expected information signals an error, as it does for
the third and fourth orders. The step of a stencil of order \\r\\ is
\\\varepsilon^{1/(r+4)}\max(1, \|x\|)\\ in the coordinate \\x\\, reduced
so that every node stays inside the domain.

The route on `scale = "link"` depends on where the parameter-scale fifth
order comes from. A family that registers its own method for this
generic (see
[`has_exact_d5expected()`](https://statmodels7.github.io/distributions7/reference/has_exact_d5expected.md))
returns the parameter scale, and the generic composes it with the lower
orders by Leibniz's rule and the Faa di Bruno formula, using the
derivatives of the inverse links up to the sixth. Every term is exact,
and the composition is correct up to rounding. For any other family the
fifth order is itself a stencil, and the stencil is taken in the linear
predictors, on the link-scale fourth derivative, which
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md)
computes exactly from the analytic orders and the derivatives of the
inverse links. Composing a stencil on the parameter scale instead was
measured against 50-digit references and rejected: where a parameter is
small the terms of the composition cancel, and the error of the stencil
is amplified up to a relative error of 6e-3, against 2e-9 or better for
the stencil in \\\eta\\. The default method reads `scale` and returns
the scale that it is given; a method of a family does not read it.

With an object as `method`, the result is the derivative of the matrix
that the method defines, computed as for the lower orders by Leibniz's
formula for
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
and with the moved density for
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md).
It reads observed derivatives up to the sixth
([`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md)),
so only `identity = "score"` is available; the Hessian identity would
need the seventh order and signals an error. On the link scale the
polynomial is evaluated with the observed derivatives on the link scale,
for the same reason: the identities and the moved density hold in any
coordinates.

## See also

[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md)
for the order below,
[`d5expected_names()`](https://statmodels7.github.io/distributions7/reference/d5expected_names.md),
[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md),
[`has_exact_d5expected()`](https://statmodels7.github.io/distributions7/reference/has_exact_d5expected.md)

## Examples

``` r
d <- poisson_distrib()
str(distrib_d5expected_hessian(d, 2, list(mu = 3), scale = "link"))
#> List of 1
#>  $ mu_mu_mu_mu_mu_mu_mu: num -3

g <- gaussian1_distrib()
r <- distrib_d5expected_hessian(g, 0, list(mu = 0, sigma = 2))
r[["sigma_sigma_sigma_sigma_sigma_sigma_sigma"]]
#> [1] 11.25
```
