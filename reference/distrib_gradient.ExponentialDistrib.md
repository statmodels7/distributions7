# Exponential Score

Computes the first derivative of the exponential log-density with
respect to the mean, one value per observation, in closed form:
\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}.\$\$
The log-density is \\-\log\mu - y/\mu\\, so the score is proportional to
the residual and its sum vanishes exactly at \\\hat\mu = \bar y\\. The
family has one parameter, so the returned list has one component.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations, non-negative.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must be strictly positive.

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

A named list of one numeric vector, `mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell\\ is the log-density of one observation and \\\mu \> 0\\ the
mean, which is also the standard deviation of this family.

## See also

[`distrib_hessian.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ExponentialDistrib.md)
for the second derivative,
[`distrib_expected_hessian.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ExponentialDistrib.md)
for the information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
y <- c(0.3, 1.1, 4.0)
th <- list(mu = 2)

# The closed form, written out.
all.equal(distrib_gradient(d, y, th)$mu, (y - 2) / 2^2)
#> [1] TRUE

# One component, the family having one parameter.
names(distrib_gradient(d, y, th))
#> [1] "mu"

# The summed score vanishes at the sample mean, which is the estimate.
set.seed(21)
z <- distrib_rng(d, 1000, list(mu = 3))
sum(distrib_gradient(d, z, list(mu = mean(z)))$mu)
#> [1] -2.287667e-14

# On the link scale the component is multiplied by h' = mu, the derivative
# of the inverse log link.
distrib_gradient(d, y, th, scale = "link")$mu / distrib_gradient(d, y, th)$mu
#> [1] 2 2 2
```
