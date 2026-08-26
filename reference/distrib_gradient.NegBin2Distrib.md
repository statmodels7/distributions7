# Negative Binomial Score, NB2

Computes the first derivatives of the negative binomial log-mass with
respect to \\\mu\\ and \\\theta\\, one value per observation. Writing
\\s = \theta + \mu\\ and \\\psi\\ for the digamma function,
\$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\theta}{s}\left(\dfrac{y}{\mu} - 1\right), \qquad \dfrac{\partial
\ell}{\partial \theta} = \psi(y+\theta) - \psi(\theta) +
\log\dfrac{\theta}{s} + \dfrac{\mu - y}{s}.\$\$ The mean component is
the score of a generalized linear model with a quadratic variance
function.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- y:

  A numeric vector of counts.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

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

A named list of two numeric vectors, `mu` and `theta`, each of length
`max(length(y), length(mu), length(theta))`.

## The dispersion at large theta

The dispersion component is **not** computed as written above. As
\\\theta\\ grows the family tends to the Poisson and this derivative
vanishes, its three terms canceling to leading order: \\\psi(y+\theta) -
\psi(\theta)\\ is \\y/\theta\\, \\\log\\\theta/s\\\\ is \\-\mu/\theta\\
and \\(\mu-y)/s\\ is \\(\mu-y)/\theta\\, and the three sum to zero, so
the value is \\O(\theta^{-2})\\ computed from terms of size
\\\theta^{-1}\\. Written directly it is wrong by about one part in a
thousand at \\\theta = 10^6\\ and **changes sign** at \\10^8\\.

The kernel performs each cancellation symbolically instead. With \\a =
\theta\\, \\b = \theta + y\\, \\c = \theta + \mu\\ and \\w =
(y-\mu)/c\\, \$\$\dfrac{\partial \ell}{\partial \theta} =
\left\\\psi(b) - \psi(a) - \log(1 + y/a)\right\\ + \left\\\log(1 + w) -
w\right\\,\$\$ each bracket evaluated by its own series where the
arguments are large. To leading order the result is \\\\y -
(y-\mu)^2\\/(2\theta^2)\\, which the example below reproduces. A fit
reaches that regime routinely: 2,000 counts drawn at a true \\\theta\\
of 100 can report an estimate of order \\10^7\\, and where such a run
stops is decided by this arithmetic.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \> 0\\ the mean and
\\\theta \> 0\\ the dispersion, with \\\operatorname{Var}(Y) = \mu +
\mu^2/\theta\\. \\\psi\\ is the digamma function, \\\psi =
(\log\Gamma)'\\.

## See also

[`distrib_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin2Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin2Distrib.md)
for their expectation and for what it does at large \\\theta\\, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 2)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out at a dispersion where the direct
# expression still has its digits.
s <- 2 + 4
all.equal(g$mu, (2 / s) * (y / 4 - 1))
#> [1] TRUE
all.equal(g$theta, digamma(y + 2) - digamma(2) + log(2 / s) + (4 - y) / s)
#> [1] TRUE

# At large theta the direct form loses them and the kernel does not. The
# value is {y - (y - mu)^2}/(2 theta^2) to leading order; at y = 3, mu = 4
# that is 1/theta^2.
cmp <- function(t) {
  c(kernel = distrib_gradient(d, 3, list(mu = 4, theta = t))$theta,
    direct = digamma(3 + t) - digamma(t) + log(t / (t + 4)) + 1 / (t + 4),
    leading = 1 / t^2)
}
rbind(`1e+04` = cmp(1e4), `1e+06` = cmp(1e6), `1e+08` = cmp(1e8))
#>             kernel        direct leading
#> 1e+04 9.999666e-09  9.999664e-09   1e-08
#> 1e+06 9.999997e-13  1.001075e-12   1e-12
#> 1e+08 1.000000e-16 -6.951360e-16   1e-16

# Summed over a fitted sample the score is at the optimizer's tolerance.
set.seed(5)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
#>           mu        theta 
#> 4.972585e-14 7.359734e-04 
```
