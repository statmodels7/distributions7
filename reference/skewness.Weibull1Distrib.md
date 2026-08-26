# Skewness of the Weibull Distribution

Closed form: \$\$\gamma_1 = \frac{g_3 - 3 g_1 g_2 + 2 g_1^3}{(g_2 -
g_1^2)^{3/2}}, \qquad g_k = \Gamma(1 + k/\sigma).\$\$ The scale cancels,
so the value is a function of the shape alone. It falls through zero at
a shape near 3.60235, so a Weibull is right-skewed below that and
left-skewed above it, which is one of the reasons the family covers more
shapes than a gamma.

## Arguments

- x:

  A `Weibull1Distrib`, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive), each a numeric vector of length 1 or `n`. The
  third moment requires a shape above \\1/3\\.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length
`max(length(theta$mu), length(theta$sigma))`. Only the shape enters the
value, so a setting that varies the scale alone repeats one number.

## Notation

\\\sigma \> 0\\ is the shape and \\g_k = \Gamma(1 + k/\sigma)\\.

## See also

[`kurtosis.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Weibull1Distrib.md),
also free of the scale;
[`variance.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
which is not;
[`weibull_gamma_factors()`](https://statmodels7.github.io/distributions7/reference/weibull_gamma_factors.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()

# At shape 1 the family is exponential, whose skewness is 2.
all.equal(skewness(d, list(mu = 1, sigma = 1)), 2)
#> [1] TRUE

# The scale does not enter a standardized moment.
skewness(d, list(mu = c(0.1, 1, 100), sigma = 2))
#> [1] 0.6311107 0.6311107 0.6311107

# The sign changes at a shape of about 3.60235.
round(skewness(d, list(mu = 1, sigma = c(2, 3.60235, 10))), 6)
#> [1]  0.631111  0.000000 -0.637637
```
