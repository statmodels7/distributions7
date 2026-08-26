# Negative Binomial Log-CDF Gradient

Closed form in the mean and an exact sum in the dispersion. The mean
component collapses like the Poisson's, \$\$\frac{\partial
F(k)}{\partial\mu} = -f(k)\\\frac{k + \theta}{\theta + \mu},\$\$ which
tends to the Poisson's \\-f(k)\\ as \\\theta\\ grows. The dispersion
component is a derivative of the incomplete beta function with respect
to its parameter, which has no elementary form, so it keeps the exact
summation of
[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md).

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- q:

  A numeric vector of quantiles. Non-integer values are floored; values
  below zero give a derivative of zero in both components.

- theta:

  A named list with components `mu` (positive) and `theta` (positive),
  each a numeric vector of length 1 or `n`. The cost of the dispersion
  component grows with the largest quantile, the sum running the support
  up to it.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `theta`, each the length
of `q` recycled against the parameters.

## Notation

\\\mu \> 0\\ is the mean, \\\theta \> 0\\ the dispersion, \\f\\ the mass
function, \\F\\ the distribution function and \\k = \lfloor q \rfloor\\.

## See also

[`discrete_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/discrete_cdf_deriv.md)
for the summation the dispersion uses;
[`distrib_grad_cdf.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PoissonDistrib.md)
for the limit;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

## Examples

``` r
d <- negbin2_distrib()
q <- c(2, 5)
th <- list(mu = 3, theta = 2)

# The mean component, written out.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th) * (q + 2) / (2 + 3))
#> [1] TRUE

# It tends to the Poisson's -f(k) as the dispersion grows.
rbind(negbin = distrib_grad_cdf(d, q, list(mu = 3, theta = 1e8),
                                log = FALSE)$mu,
      poisson = distrib_grad_cdf(poisson_distrib(), q, list(mu = 3),
                                 log = FALSE)$mu)
#>               [,1]       [,2]
#> negbin  -0.2240418 -0.1008188
#> poisson -0.2240418 -0.1008188
```
