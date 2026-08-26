# Chi-Squared Observed Hessian

Computes the second derivative of the chi-squared log-density with
respect to \\\mu\\, in closed form: \$\$\dfrac{\partial^2 \ell}{\partial
\mu^2} = -\dfrac{\psi'(\mu/2)}{4},\$\$ with \\\psi'\\ the trigamma
function. **It does not involve the response.** The family is a
one-parameter exponential family in \\\log y\\, so the data reach the
log-density only through a term linear in \\\mu\\ and the second
derivative kills it. The value is constant within a parameter setting,
is recycled to the length of `y`, and equals
[`distrib_expected_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ChisqDistrib.md)
exactly.

The coincidence is on the **parameter** scale. On the link scale the two
differ, the second-order chain rule adding a term
\\h''(\eta)\\\partial\ell/\partial\mu\\ that the expected version drops
and a sample does not; the example below measures it.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

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

A named list with one numeric vector, `mu_mu`, of length
`max(length(y), length(mu))` and constant within itself when the
parameter is.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.
\\\psi'\\ is the trigamma function and \\h(\eta)\\ the inverse link.

## See also

[`distrib_gradient.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ChisqDistrib.md)
for the score, which is the one quantity here that reads the data;
[`distrib_expected_hessian.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ChisqDistrib.md),
which returns the same number on the parameter scale;
[`distrib_deriv3.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ChisqDistrib.md)
for the order above; and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
y <- c(1, 4, 9)
th <- list(mu = 4)

# A constant across the observations, written out with the trigamma.
distrib_hessian(d, y, th)
#> $mu_mu
#> [1] -0.1612335 -0.1612335 -0.1612335
#> 
-trigamma(2) / 4
#> [1] -0.1612335

# Free of the response, so it equals its own expectation to the bit.
identical(distrib_hessian(d, y, th), distrib_expected_hessian(d, y, th))
#> [1] TRUE

# On the link scale they differ by h''(eta) times the score, which is what
# makes Fisher scoring and Newton's method take different steps.
obs <- distrib_hessian(d, y, th, scale = "link")$mu_mu
exp_ <- distrib_expected_hessian(d, y, th, scale = "link")$mu_mu
rbind(observed = obs, expected = exp_,
      difference = obs - exp_,
      h2_times_score = 4 * distrib_gradient(d, y, th)$mu)
#>                     [,1]       [,2]       [,3]
#> observed       -4.811599 -2.0390106 -0.4171501
#> expected       -2.579736 -2.5797363 -2.5797363
#> difference     -2.231863  0.5407257  2.1625861
#> h2_times_score -2.231863  0.5407257  2.1625861
```
