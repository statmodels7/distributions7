# Gradient of the Log Distribution Function

Computes the first derivatives, with respect to the parameters, of
\\\log F(q;\theta)\\ — or of \\\log(1 - F(q;\theta))\\ when
`lower.tail = FALSE`.

These are what a **censored** observation contributes to the score. An
observation known only to be at most \\q\\ contributes \\\log F(q)\\,
one known only to exceed \\q\\ contributes \\\log(1-F(q))\\, and an
interval censored one contributes \\\log(F(b) - F(a))\\, which is
assembled from the unlogged derivatives (`log = FALSE`) at the two
endpoints. They are also what the delta method needs for the standard
error of a quantile residual.

## Usage

``` r
distrib_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)
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

A named list with one numeric vector per parameter.

## Details

The mathematics is one exchange of derivative and integral. Since the
region of integration does not depend on \\\theta\\, \$\$\frac{\partial
F(q;\theta)}{\partial\theta_i} = \int\_{-\infty}^{q}\frac{\partial
f}{\partial\theta_i} = \int\_{-\infty}^{q} f\\\ell^{(i)} =
F(q)\\\mathbb{E}\\\left\[\ell^{(i)} \mid Y \leq q\right\],\$\$ so the
gradient of the log distribution function is a *partial mean of the
score*. For a discrete distribution the integral is a finite sum and the
identity is exact, which is how the default method computes it there;
for a continuous one the default differentiates
[`distrib_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
numerically, and distributions with a closed form register it directly.

## See also

[`distrib_hess_cdf`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
[`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)

## Examples

``` r
d <- gaussian_distrib()
theta <- list(mu = 0, sigma = 1)

# what a right-censored observation at q = 1 contributes to the score
distrib_grad_cdf(d, 1, theta, lower.tail = FALSE)
#> $mu
#> [1] 1.525135
#> 
#> $sigma
#> [1] 1.525135
#> 
```
