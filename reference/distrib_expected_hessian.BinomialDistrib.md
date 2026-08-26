# Binomial Expected Hessian

Returns the expectation of the observed second derivative under the
model, in closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{n}{\mu(1-\mu)},\$\$ which follows from \\\mathbb{E}\[Y\] =
n\mu\\ in \\-y/\mu^2 - (n-y)/(1-\mu)^2\\. The Fisher information is
\\n/(\mu(1-\mu))\\, so it grows in proportion to the number of trials
and is smallest at \\\mu = 1/2\\.

On the **link** scale with the default logit the value is
\\-n\mu(1-\mu)\\, and so is the observed Hessian: the logit is the
canonical link and the two coincide exactly, so iteratively reweighted
least squares is Fisher scoring and Newton's method at once.

Because the value does not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials.

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(0, 1)\\.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available. Accepted so that the signature
  matches the generic's.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu`, of length
`max(length(y), length(mu), length(distrib@size))`, constant at
\\-n/(\mu(1-\mu))\\ where `size` is.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\mu^2\]\\, the expectation of the
**observed information** under the model. The binomial is a regular
family, so the second Bartlett identity holds and this equals the
variance of the score.

## See also

[`distrib_hessian.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BinomialDistrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)
th <- list(mu = 0.3)

# A single number, -n divided by the Bernoulli variance.
unique(distrib_expected_hessian(d, c(0, 4, 10), th)$mu_mu)
#> [1] -47.61905
-10 / (0.3 * 0.7)
#> [1] -47.61905

# It grows in proportion to the number of trials.
vapply(c(1, 10, 100), function(n)
  -distrib_expected_hessian(binomial_distrib(size = n), 0, th)$mu_mu,
  numeric(1))
#> [1]   4.761905  47.619048 476.190476

# On the canonical logit link the observed and expected values agree
# exactly, at every observation.
rbind(observed = distrib_hessian(d, c(0, 4, 10), th, scale = "link")$mu_mu,
      expected = distrib_expected_hessian(d, c(0, 4, 10), th,
                                          scale = "link")$mu_mu)
#>          [,1] [,2] [,3]
#> observed -2.1 -2.1 -2.1
#> expected -2.1 -2.1 -2.1
```
