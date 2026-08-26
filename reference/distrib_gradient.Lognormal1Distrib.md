# Lognormal Score

Computes the first derivatives of the lognormal log-density with respect
to \\\mu\\ and \\\sigma^2\\, one value per observation, in closed form.
With \\r = \log y - \mu\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{r}{\sigma^2}, \qquad \dfrac{\partial \ell}{\partial \sigma^2} =
\dfrac{r^2 - \sigma^2}{2\sigma^4}.\$\$

**These are exactly the Gaussian's, read at \\\log y\\.** The lognormal
log-density is the Gaussian's at \\\log y\\ minus \\\log y\\, and that
last term carries no parameter, so it disappears from every derivative
in \\\mu\\ and \\\sigma^2\\. The same holds at all four orders and for
the expected values; what differs from
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md)
is only the derivatives in the response, where the Jacobian does enter.

The arithmetic runs in a compiled kernel decomposed over the elements of
the output, so the result does not depend on the thread count.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of two numeric vectors, `mu` and `sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and
\\\sigma^2 \> 0\\ the variance **of \\\log Y\\**, not of \\Y\\. \\r =
\log y - \mu\\ is the residual on the log scale.

## See also

[`distrib_hessian.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Lognormal1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Lognormal1Distrib.md)
for their expectation,
[`distrib_gradient.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian2Distrib.md),
which returns the same numbers at \\\log y\\, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out on the log scale.
r <- log(y) - 0.5
all.equal(g$mu, r / 0.36)
#> [1] TRUE
all.equal(g$sigma2, (r^2 - 0.36) / (2 * 0.36^2))
#> [1] TRUE

# Identical to the Gaussian's score at log y, component for component.
all.equal(g, distrib_gradient(gaussian2_distrib(), log(y), th))
#> [1] TRUE

# The summed score vanishes at the closed-form estimates, which are the
# sample moments of the logarithm.
set.seed(6)
z <- distrib_rng(d, 2000, th)
mle <- list(mu = mean(log(z)),
            sigma2 = mean((log(z) - mean(log(z)))^2))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>           mu       sigma2 
#> 1.814434e-13 6.810524e-15 
```
