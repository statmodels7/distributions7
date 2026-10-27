# Chi-Squared Third-Order Derivative

Computes the single third derivative of the chi-squared log-density with
respect to \\\mu\\, in closed form: \$\$\ell^{(\mu\mu\mu)} =
-\dfrac{\psi''(\mu/2)}{8},\$\$ with \\\psi''\\ the second derivative of
the digamma function. It is one case of the general pattern \\\ell^{(k)}
= -\psi^{(k-2)}(\mu/2)/2^{k}\\ for \\k \ge 2\\, which holds because the
family is a one-parameter exponential family in \\\log y\\ and the
response leaves the derivatives after the first order.

On the parameter scale the value is free of the response, so `expected`
selects nothing: the same kernel runs either way and the two results are
identical to the bit. With `scale = "link"` the generic composes the
observed and the expected derivative differently, and the two results
differ.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`, recycled if of length 1. It must be strictly
  positive.

- expected:

  Logical of length 1. On the parameter scale it has no effect, because
  the observed and expected third derivatives are the same number; on
  the link scale it selects the expected composition. Defaults to
  `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count is. Defaults to `1L`.

## Value

A named list with one numeric vector, `mu_mu_mu`, of length
`max(length(y), length(mu))` and constant within itself when the
parameter is.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. \\\psi\\ is the digamma function and
\\\psi^{(m)}\\ its \\m\\th derivative.

## See also

[`distrib_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ChisqDistrib.md)
for the order below and
[`distrib_deriv4.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ChisqDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
y <- c(1, 4, 9)
th <- list(mu = 4)

# The value is a constant and is written out with the second derivative
# of the digamma function.
unique(distrib_deriv3(d, y, th)$mu_mu_mu)
#> [1] 0.05051423
-psigamma(2, deriv = 2) / 8
#> [1] 0.05051423

# Free of the response, so on the parameter scale the expectation is the
# same.
identical(distrib_deriv3(d, y, th), distrib_deriv3(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the Hessian reproduces it.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 4 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 4 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
          tolerance = 1e-6)
#> [1] TRUE
```
