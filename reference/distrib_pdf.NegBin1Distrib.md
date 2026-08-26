# NB1 Probability Mass Function

Computes the negative binomial mass at size \\r = \mu/\theta\\ and
success probability \\1/(1+\theta)\\, \$\$P(Y = y; \mu, \theta) =
\dfrac{\Gamma(y+r)}{\Gamma(r)\\y!} \left(\dfrac{1}{1+\theta}\right)^{r}
\left(\dfrac{\theta}{1+\theta}\right)^{y}, \qquad y = 0, 1, 2,
\ldots,\$\$ which is the pairing that makes the variance
\\\mu(1+\theta)\\. The log-mass is formed in a compiled kernel and
exponentiated when `log` is `FALSE`, so `log = TRUE` is the accurate
route far into the tail.

The mean sits **inside** the gamma functions here, through the size. In
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
the size is \\\theta\\ and the mean stays outside them, which is the
concrete difference between the two families.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- y:

  A numeric vector of counts. A negative value gives a mass of 0.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned. Defaults to
  `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities, of length
`max(length(y), length(mu), length(theta))`, one value per observation.

## See also

[`distrib_cdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin1Distrib.md)
for the distribution function,
[`distrib_gradient.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin1Distrib.md)
for the derivatives of the log-mass,
[`nb1_size()`](https://statmodels7.github.io/distributions7/reference/nb1_size.md)
and
[`nb1_prob()`](https://statmodels7.github.io/distributions7/reference/nb1_prob.md)
for the pairing this uses, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 4)

# The mass is stats::dnbinom at size mu/theta and prob 1/(1 + theta).
all.equal(distrib_pdf(d, y, th), dnbinom(y, size = 1, prob = 1 / 5))
#> [1] TRUE

# At mu = theta the size is 1 and the family is the geometric.
all.equal(distrib_pdf(d, 0:4, th), dgeom(0:4, prob = 1 / 5))
#> [1] TRUE

# As theta goes to zero the mass tends to the Poisson's.
rbind(nb1 = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e-6)),
      poisson = dpois(0:4, 3))
#>               [,1]      [,2]      [,3]      [,4]      [,5]
#> nb1     0.04978714 0.1493613 0.2240418 0.2240417 0.1680313
#> poisson 0.04978707 0.1493612 0.2240418 0.2240418 0.1680314

# The mass sums to one over the support.
sum(distrib_pdf(d, 0:600, th))
#> [1] 1
```
