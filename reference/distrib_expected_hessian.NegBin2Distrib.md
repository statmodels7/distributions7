# Negative Binomial Expected Hessian, NB2

Returns the expectation of the observed Hessian under the model. With
\\s = \theta + \mu\\ and \\\psi_1\\ the trigamma function,
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{\theta}{\mu s},
\qquad \mathbb{E}\left\[\ell^{(\mu\theta)}\right\] = 0, \qquad
\mathbb{E}\left\[\ell^{(\theta\theta)}\right\] =
\mathbb{E}\[\psi_1(Y+\theta)\] - \psi_1(\theta) + \dfrac{\mu}{\theta
s}.\$\$ The mean entry and the mixed one are closed form, the latter
vanishing because the observed mixed entry is \\(y-\mu)/s^2\\ and
\\\mathbb{E}\[Y\] = \mu\\. So the mean and the dispersion are
orthogonal.

\\\mathbb{E}\[\psi_1(Y+\theta)\]\\ has no closed form. It is summed over
the support through the mass recurrence, carried in log scale until the
terms are representable, and stopped when the accumulated mass reaches
\\1 - 10^{-12}\\. The sum is exact to that mass, so this entry is a
truncated exact sum and not a quadrature or a simulation; `approx` and
`nsim` are ignored, and `y` is read only for its length.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- y:

  A numeric vector of counts. Only its length is used.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored. The mean and mixed entries are closed form and the dispersion
  entry is a truncated exact sum. Accepted so that the signature matches
  the generic's.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `theta_theta` and
`mu_theta`, in that order, each of length
`max(length(y), length(mu), length(theta))` and constant within itself
when the parameters are.

## A caveat at large theta

The dispersion entry is a difference of two quantities that agree to
leading order as \\\theta\\ grows, and **it is not rewritten to remove
that cancellation**, unlike the score and the observed Hessian. Its
leading order needs one term more of the observed Hessian than those two
do, the \\\theta^{-3}\\ term vanishing under expectation, so it is a
derivation of its own. Measured at \\\mu = 4\\: it reads \\-7.3\times
10^{-8}\\ at \\\theta = 10^2\\, \\-3.6\times 10^{-16}\\ at \\10^4\\ and
\\+1.7\times 10^{-16}\\ at \\10^6\\. The last is **positive**, which an
expected second derivative cannot be, so a Fisher scoring step taken
there is not reliable. The example below shows the sign turning.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
negative binomial is a regular family, so the second Bartlett identity
holds and this equals the variance of the score.

## See also

[`distrib_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin2Distrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()
th <- list(mu = 4, theta = 2)

# The three entries, one value per observation.
lapply(distrib_expected_hessian(d, c(0, 2, 6), th), unique)
#> $mu_mu
#> [1] -0.08333333
#> 
#> $theta_theta
#> [1] -0.0492447
#> 
#> $mu_theta
#> [1] 0
#> 

# The mean entry is closed form and the mixed one is exactly zero.
c(closed = -2 / (4 * (2 + 4)),
  mixed = distrib_expected_hessian(d, 0, th)$mu_theta)
#>      closed       mixed 
#> -0.08333333  0.00000000 

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>                mu_mu theta_theta    mu_theta
#> observed -0.08474479 -0.04931509 0.001129167
#> expected -0.08333333 -0.04924470 0.000000000

# The dispersion entry cancels at large theta and turns positive, which an
# expected second derivative cannot be.
vapply(c(1e2, 1e4, 1e6),
       function(t) distrib_expected_hessian(d, 0,
                     list(mu = 4, theta = t))$theta_theta,
       numeric(1))
#> [1] -7.326902e-08 -3.581440e-16  1.734243e-16
```
