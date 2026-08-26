# Chi-Squared Score

Computes the derivative of the chi-squared log-density with respect to
\\\mu\\, one value per observation, in closed form: \$\$\dfrac{\partial
\ell}{\partial \mu} = \dfrac{\log y - \log 2 - \psi(\mu/2)}{2},\$\$ with
\\\psi\\ the digamma function. The family is a one-parameter exponential
family in \\\log y\\, so this is the sufficient statistic minus its
expectation, and **it is the only order that involves the response at
all**: every derivative from the second up is a polygamma function of
\\\mu/2\\ and nothing else.

The score has mean zero because \\\mathbb{E}\[\log Y\] = \psi(\mu/2) +
\log 2\\.

With `scale = "link"` the generic applies the chain rule for the link
the family carries before returning. This method always returns the
parameter scale.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of strictly positive observations. A value at zero
  makes the logarithm infinite and the score non-finite.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`, recycled if of length 1. It must be strictly
  positive.

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

A named list with one numeric vector, `mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell\\ is the log-density of one observation and \\\mu \> 0\\ the
mean, which is also the degrees of freedom. \\\psi\\ is the digamma
function, \\\psi = (\log\Gamma)'\\.

## See also

[`distrib_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ChisqDistrib.md)
for the second derivative,
[`distrib_expected_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ChisqDistrib.md),
which returns the same number,
[`distrib_grad_y.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.ChisqDistrib.md)
for the derivative in the response, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
y <- c(1, 4, 9)
th <- list(mu = 4)
g <- distrib_gradient(d, y, th)

# The closed form, written out with the digamma function.
all.equal(g$mu, (log(y) - log(2) - digamma(2)) / 2)
#> [1] TRUE

# It is a sufficient statistic minus its expectation, so the sample mean of
# log y matches psi(mu/2) + log 2.
set.seed(9)
z <- distrib_rng(d, 2e5, th)
c(sample = mean(log(z)), theory = digamma(2) + log(2))
#>   sample   theory 
#> 1.116566 1.115932 

# The score vanishes where log y equals that expectation.
distrib_gradient(d, exp(digamma(2) + log(2)), th)$mu
#> [1] -5.551115e-17

# Summed over a fitted sample it is at the optimizer's tolerance.
set.seed(7)
zz <- distrib_rng(d, 2000, th)
sum(distrib_gradient(d, zz, as.list(coef(fit_distrib(d, zz))))$mu)
#> [1] 3.224567e-10
```
