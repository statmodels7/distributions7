# Negative Binomial Probability Mass Function, NB2

Computes the negative binomial mass \$\$P(Y = y; \mu, \theta) =
\dfrac{\Gamma(y+\theta)}{y!\\\Gamma(\theta)}
\left(\dfrac{\theta}{\theta+\mu}\right)^{\theta}
\left(\dfrac{\mu}{\theta+\mu}\right)^{y}, \qquad y = 0, 1, 2,
\ldots,\$\$ by calling
[`stats::dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) at
`size = theta` and `mu = mu`, so the accuracy is R's own. With
`log = TRUE` the logarithm is formed inside
[`dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) and stays finite
far into the tail.

The mass is that of a Poisson whose rate is itself gamma distributed
with mean \\\mu\\ and shape \\\theta\\, which is where the
overdispersion comes from. As \\\theta\\ grows the gamma concentrates
and the mass tends to the Poisson's; at \\\theta = 1\\ it is the
geometric.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- y:

  A numeric vector of counts. A non-integer value gives 0 with a warning
  from [`stats::dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html),
  and a negative value gives 0.

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

[`distrib_cdf.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.NegBin2Distrib.md)
for the distribution function,
[`distrib_gradient.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin2Distrib.md)
for the derivatives of the log-mass,
[`distrib_pdf.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.NegBin1Distrib.md)
for the NB1 parametrization, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 2)

# The method is stats::dnbinom at size = theta, mu = mu.
all.equal(distrib_pdf(d, y, th), dnbinom(y, size = 2, mu = 4))
#> [1] TRUE

# At theta = 1 the family is the geometric.
all.equal(distrib_pdf(d, y, list(mu = 4, theta = 1)),
          distrib_pdf(geometric_distrib(), y, list(mu = 4)))
#> [1] TRUE

# At a large theta it is the Poisson, to six figures.
rbind(negbin = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e6)),
      poisson = dpois(0:4, 3))
#>               [,1]      [,2]      [,3]      [,4]      [,5]
#> negbin  0.04978729 0.1493614 0.2240417 0.2240415 0.1680311
#> poisson 0.04978707 0.1493612 0.2240418 0.2240418 0.1680314

# The mass sums to one over the support.
sum(distrib_pdf(d, 0:400, th))
#> [1] 1
```
