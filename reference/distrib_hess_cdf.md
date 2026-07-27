# Second Derivatives of the Log Distribution Function

Computes the second derivatives, with respect to the parameters, of
\\\log F(q;\theta)\\ — or of \\\log(1 - F(q;\theta))\\ when
`lower.tail = FALSE`. Together with
[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
these give the observed information of a censored observation.

## Usage

``` r
distrib_hess_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- q:

  A numeric vector of quantiles.

- theta:

  A named list (or named numeric vector) of distribution parameters.
  Each parameter must have length 1 or `length(q)`.

- lower.tail:

  Logical; if `TRUE` (default), derivatives of \\\log F(q)\\, otherwise
  of \\\log(1 - F(q))\\.

- log:

  Logical; if `TRUE` (default), derivatives of the *log* tail
  probability. With `FALSE` the derivatives of the probability itself
  are returned, which is what interval censoring and the truncation
  constant are built from.

- ...:

  Additional arguments passed to the specific method.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)`(distrib@params)`.

## Details

By the same exchange as in
[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
and using \\\partial\_{ij} f / f = \ell^{(ij)} + \ell^{(i)}\ell^{(j)}\\,
\$\$\frac{\partial^{2} F(q)}{\partial\theta_i\partial\theta_j} =
F(q)\\\mathbb{E}\\\left\[\ell^{(ij)} + \ell^{(i)}\ell^{(j)} \mid Y \leq
q\right\],\$\$ and the log scale follows from \\\partial\_{ij}\log P =
\partial\_{ij}P/P - (\partial_i P/P)(\partial_j P/P)\\.

## See also

[`distrib_grad_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
