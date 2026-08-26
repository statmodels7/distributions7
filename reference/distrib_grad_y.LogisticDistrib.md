# Logistic First Derivative in the Response

Computes the first derivative of the logistic log-density with respect
to the response, \$\$\dfrac{\partial \ell}{\partial y} =
-\dfrac{1}{\sigma}\tanh\left(\dfrac{z}{2}\right), \qquad z =
\dfrac{y-\mu}{\sigma},\$\$ in closed form. The logistic is a location
family in \\\mu\\, so the response enters the log-density only through
\\z\\ and this derivative is the negative of the score in \\\mu\\. It is
bounded by \\1/\sigma\\.

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

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu), length(sigma))`,
one value per observation.

## See also

[`distrib_hess_y.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.LogisticDistrib.md)
for the second derivative in the response,
[`distrib_gradient.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LogisticDistrib.md)
for the score in the parameters, and
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)

all.equal(distrib_grad_y(d, y, th), -tanh((y - 0.4) / (2 * 1.5)) / 1.5)
#> [1] TRUE

# A location family: the derivative in the response is minus the score in
# the location.
all.equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
#> [1] TRUE

# Bounded by 1/sigma however far out the observation is.
max(abs(distrib_grad_y(d, seq(-1e3, 1e3, length.out = 1e4), th)))
#> [1] 0.6666667
1 / 1.5
#> [1] 0.6666667
```
