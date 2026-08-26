# Gaussian Score in Mean and Variance

Computes the first derivatives of the Gaussian log-density with respect
to \\\mu\\ and \\v = \sigma^2\\, one value per observation, in closed
form. With \\r = y - \mu\\, \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{r}{v}, \qquad \dfrac{\partial \ell}{\partial v} = \dfrac{r^2 -
v}{2v^2}.\$\$ The arithmetic runs in a compiled kernel decomposed over
the elements of the output, so the result does not depend on the thread
count.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning, giving \\\partial \ell / \partial
\eta_j = h_j'(\eta_j)\\ \partial \ell / \partial \theta_j\\. This method
always returns the parameter scale; the transformation happens in the
generic.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- y:

  A numeric vector of observations.

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

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean and \\v
= \sigma^2 \> 0\\ the variance. \\\eta_j\\ is the coordinate of
parameter \\j\\ on the unconstrained scale of its link, and \\h_j' =
\partial \theta_j / \partial \eta_j\\ the chain-rule factor onto it.

## See also

[`distrib_hessian.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian2Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian2Distrib.md)
for their expectation,
[`distrib_gradient.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian1Distrib.md)
for the same score in the standard deviation, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gaussian2_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 1, sigma2 = 4)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
r <- y - 1
all.equal(g$mu, r / 4)
#> [1] TRUE
all.equal(g$sigma2, (r^2 - 4) / (2 * 4^2))
#> [1] TRUE

# The summed score vanishes at the maximum likelihood estimate, where the
# variance carries the divisor n.
set.seed(1)
z <- rnorm(200, mean = 3, sd = 2)
mle <- list(mu = mean(z), sigma2 = mean((z - mean(z))^2))
vapply(distrib_gradient(d, z, mle), sum, numeric(1))
#>            mu        sigma2 
#> -5.859029e-15 -1.425943e-15 

# Both parametrizations ride a log link on the spread, and log(sigma^2) is
# twice log(sigma), so the link-scale score here is half gaussian1's.
distrib_gradient(d, y, th, scale = "link")$sigma2 /
  distrib_gradient(gaussian1_distrib(), y, list(mu = 1, sigma = 2),
                   scale = "link")$sigma
#> [1] 0.5 0.5 0.5
```
