# Chi-Squared Second Derivative in the Response

Computes the second derivative of the chi-squared log-density with
respect to the response, in closed form: \$\$\dfrac{\partial^2
\ell}{\partial y^2} = -\dfrac{\mu/2 - 1}{y^2}.\$\$ The exponential term
drops out, the log-density being linear in \\y\\ apart from the
\\(\mu/2 - 1)\log y\\ term. The sign follows \\\mu\\: the log-density is
concave in the response for \\\mu \> 2\\, exactly flat at \\\mu = 2\\,
where the family is the exponential with mean 2, and convex for \\\mu \<
2\\.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of strictly positive observations. At `y = 0` the
  value is infinite unless \\\mu\\ is exactly 2.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`, recycled if of length 1. It must be strictly
  positive.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of length `max(length(y), length(mu))`, one value per
observation.

## See also

[`distrib_grad_y.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ChisqDistrib.md)
for the first derivative in the response,
[`distrib_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ChisqDistrib.md)
for the curvature in the parameter, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
y <- c(1, 4, 9)
th <- list(mu = 4)

all.equal(distrib_hess_y(d, y, th), -(4 / 2 - 1) / y^2)
#> [1] TRUE

# Concave above mu = 2, flat at it, convex below.
vapply(c(4, 2, 1), function(m) distrib_hess_y(d, 4, list(mu = m)),
       numeric(1))
#> [1] -0.06250  0.00000  0.03125

# Negative everywhere at four degrees of freedom.
all(distrib_hess_y(d, y, th) < 0)
#> [1] TRUE
```
