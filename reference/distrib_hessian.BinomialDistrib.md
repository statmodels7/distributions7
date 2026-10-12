# Binomial Observed Hessian

Computes the second derivative of the binomial log-mass with respect to
the probability, one value per observation, in closed form:
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} -
\dfrac{n-y}{(1-\mu)^2}.\$\$ It is negative for every admissible \\\mu\\
and depends on the data through \\y\\ alone. Each term is formed by
dividing its coefficient by the base once per power, and a term with a
zero coefficient contributes zero, so at \\y = 0\\ the value is
\\-n/(1-\mu)^2\\ however small \\\mu\\ is.

On the **link** scale with the default logit the chain rule gives
\\\partial^2\ell/\partial\eta^2 = -n\mu(1-\mu)\\, which carries no data;
this is the defining property of a canonical link, and in exact
arithmetic it makes the observed and the expected information the same
matrix there (see
[`distrib_expected_hessian.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BinomialDistrib.md)).
The generic forms the link-scale value from the derivatives on the
parameter scale, and for an observation on the far side of the
probability (\\y = 0\\ with \\\mu\\ near 1, or \\y = n\\ with \\\mu\\
near 0) its two terms cancel to a value much smaller than either. The
relative error of the observed value then grows with \\\|\eta\|\\ (about
\\10^{-3}\\ at \\\eta = 27.6\\), while the expected value is computed
directly and stays exact.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials.

- y:

  A numeric vector of counts of successes.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(0, 1)\\.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu`, of length
`max(length(y), length(mu), length(distrib@size))`. Each of `y`, `mu`
and `size` must have length 1 or that length; any other length signals
an error.

## See also

[`distrib_gradient.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BinomialDistrib.md)
for the score,
[`distrib_expected_hessian.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BinomialDistrib.md)
for the expectation of this quantity, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)
y <- c(0, 4, 10)
th <- list(mu = 0.3)

# The closed form, written out.
all.equal(distrib_hessian(d, y, th)$mu_mu,
          -y / 0.3^2 - (10 - y) / 0.7^2)
#> [1] TRUE

# On the canonical logit link the curvature carries no data: one value,
# repeated, equal to -n mu (1-mu).
distrib_hessian(d, y, th, scale = "link")$mu_mu
#> [1] -2.1 -2.1 -2.1
-10 * 0.3 * 0.7
#> [1] -2.1

# A central difference of the score reproduces the parameter-scale value.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 0.3 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 0.3 - eps))$mu
all.equal((up - dn) / (2 * eps), distrib_hessian(d, y, th)$mu_mu,
          tolerance = 1e-6)
#> [1] TRUE
```
