# NB1 Score

Computes the first derivatives of the NB1 log-mass with respect to
\\\mu\\ and \\\theta\\, one value per observation, in closed form. Both
come from the chain rule through the size \\r = \mu/\theta\\. Writing
\\P = \psi(y+r) - \psi(r) - \log(1+\theta)\\ with \\\psi\\ the digamma
function, \$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{P}{\theta},
\qquad \dfrac{\partial \ell}{\partial \theta} =
-\dfrac{\mu}{\theta^2}P - \dfrac{r}{1+\theta} + \dfrac{y}{\theta} -
\dfrac{y}{1+\theta}.\$\$ The arithmetic runs in a compiled kernel
decomposed over the elements of the output, so the result does not
depend on the thread count.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

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

## The dispersion at small theta

As \\\theta\\ goes to zero the family tends to the Poisson, so the
dispersion component tends to a finite limit while its individual terms
run away: \\r = \mu/\theta\\ grows without bound and the chain rule
divides by \\\theta^2\\. The kernel computes the digamma difference in a
form that performs the cancellation symbolically, and the value
converges onto \$\$\lim\_{\theta \to 0} \dfrac{\partial \ell}{\partial
\theta} = \dfrac{(y-\mu)^2 - y}{2\mu},\$\$ which the example below
checks at three settings. What survives is about five significant
figures at \\\theta = 10^{-8}\\ and no more, the powers of \\r\\ in the
assembly carrying a cancellation of their own that is not removed.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. This method always returns the
parameter scale.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \> 0\\ the mean and
\\\theta \> 0\\ the dispersion, with \\\operatorname{Var}(Y) =
\mu(1+\theta)\\. \\r = \mu/\theta\\ is the size and \\\psi\\ the digamma
function.

## See also

[`distrib_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin1Distrib.md)
for the second derivatives,
[`distrib_expected_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin1Distrib.md)
for their expectation and for what it does at small \\\theta\\,
[`distrib_gradient.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin2Distrib.md)
for the quadratic-variance family, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 4)
g <- distrib_gradient(d, y, th)

# The two closed forms, written out.
r <- 4 / 4
P <- digamma(y + r) - digamma(r) - log(1 + 4)
all.equal(g$mu, P / 4)
#> [1] TRUE
all.equal(g$theta, -(4 / 4^2) * P - r / (1 + 4) + y / 4 - y / (1 + 4))
#> [1] TRUE

# As theta goes to zero the dispersion component converges onto
# {(y - mu)^2 - y}/(2 mu).
lim <- function(y, mu) ((y - mu)^2 - y) / (2 * mu)
cmp <- function(y, mu) {
  c(vapply(c(1e-4, 1e-6, 1e-8),
           function(t) distrib_gradient(d, y, list(mu = mu, theta = t))$theta,
           numeric(1)),
    limit = lim(y, mu))
}
rbind(`y=3,mu=4` = cmp(3, 4), `y=0,mu=3` = cmp(0, 3), `y=7,mu=2` = cmp(7, 2))
#>                                     limit
#> y=3,mu=4 -0.2499979 -0.250000 -0.25 -0.25
#> y=0,mu=3  1.4998000  1.499998  1.50  1.50
#> y=7,mu=2  4.4982922  4.499983  4.50  4.50

# Summed over a fitted sample the score is at the optimizer's tolerance.
set.seed(5)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
vapply(distrib_gradient(d, z, as.list(coef(fit))), sum, numeric(1))
#>            mu         theta 
#>  1.856601e-05 -1.845840e-05 
```
