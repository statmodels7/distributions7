# Logistic Third-Order Derivatives

Computes the four distinct third derivatives of the logistic log-density
with respect to \\\mu\\ and \\\sigma\\. The observed values are closed
form. Writing \\z = (y-\mu)/\sigma\\, \\t = 1/(1+e^{-z})\\ and \\u =
1-t\\, the log-density is \\\ell = -\log\sigma + g(z)\\ with \\g(z) =
-z - 2\log(1+e^{-z})\\, whose derivatives in \\z\\ are \$\$g_1 = 1-2t,
\quad g_2 = -2tu, \quad g_3 = -2tu(1-2t), \quad g_4 = -2tu(1-6tu).\$\$
Both parameters enter only through \\z\\, so each component is a
polynomial in \\z\\ with the \\g_j\\ as coefficients:
\$\$\dfrac{\partial^3 \ell}{\partial \mu^3} = -\dfrac{g_3}{\sigma^3},
\qquad \dfrac{\partial^3 \ell}{\partial \mu^2 \partial \sigma} =
-\dfrac{2g_2 + z g_3}{\sigma^3},\$\$ \$\$\dfrac{\partial^3
\ell}{\partial \mu \partial \sigma^2} = -\dfrac{2g_1 + 4z g_2 + z^2
g_3}{\sigma^3}, \qquad \dfrac{\partial^3 \ell}{\partial \sigma^3} =
-\dfrac{2 + 6z g_1 + 6z^2 g_2 + z^3 g_3}{\sigma^3}.\$\$

With `expected = TRUE` the values are **numerical**, unlike the orders
below them. Two of the nine expectations at this order and the next
require \\\int w^k \mathrm{sech}^4 w \tanh^2 w \\ dw\\, which has no
elementary form, so the method routes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
and `approx` and `nsim` are read. The default there is `"integrate"`,
one quadrature per component.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned, by the numerical route described above. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  One of `"integrate"` (the default), `"bartlett"`, `"mc"` or `"opg"`,
  matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read
  **only** when `expected = TRUE`, where it selects how the expectation
  is approximated.

- nsim:

  A single positive integer, the sample size used when `approx = "mc"`.
  Read only when `expected = TRUE`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use for the
  observed values. Defaults to `1L`.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell^{(i j k)}\\ is the third derivative of the log-density with
respect to parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts
name derivatives; a subscript on \\\ell\\ never does.

## See also

[`distrib_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md)
for the order below and
[`distrib_deriv4.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LogisticDistrib.md)
for the order above;
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
for the numerical route the expected values take;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# A central difference of the Hessian reproduces the observed component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#> [1] TRUE

# The expected values are numerical here: the pure-mu component is zero by
# symmetry and comes back at the quadrature's own accuracy, not at zero.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 4.336809e-19
#> 
#> $mu_mu_sigma
#> [1] 0.1481481
#> 
#> $mu_sigma_sigma
#> [1] 5.551115e-17
#> 
#> $sigma_sigma_sigma
#> [1] 1.75846
#> 
```
