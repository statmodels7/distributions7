# Pseudo-Huber Expected Hessian

Returns the expectation of the observed Hessian under the model. **There
is no closed form**, so the four components that do not vanish are
obtained by the strategy `approx` names, normally a numerical
integration of the observed Hessian against the density through
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).
The two components containing \\\mu\\ an odd number of times are then
**replaced by exact zeros**: the law is symmetric about \\\mu\\, so
\\\mathbb{E}\[r\] = \mathbb{E}\[r^3\] = 0\\ and the \\\mu\sigma\\ and
\\\mu\nu\\ entries vanish. The location is therefore orthogonal to both
other parameters, and \\\hat\mu\\ is asymptotically independent of them.

The method **improves** the approximation rather than replacing it,
which is why
[`expected_hessian_exact.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.PseudoHuberDistrib.md)
answers `FALSE`. Reading the method's owning class would say the family
writes its information out; it does not, and the cost says so: measured
at 100 observations this takes about 11 seconds, where the families that
do write it out answer in a median of 0.183 milliseconds.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- y:

  A numeric vector of observations. Its length sets the length of each
  returned component.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. `sigma` and `nu` must be
  strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  One of `"bartlett"` (the default), `"integrate"`, `"mc"` or `"opg"`,
  the strategy
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  uses. **Read here**, unlike on the families that write their
  information out.

- nsim:

  A single positive integer, the sample size when `approx = "mc"`.
  Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of six numeric vectors, `mu_mu`, `sigma_sigma`, `nu_nu`,
`mu_sigma`, `mu_nu` and `sigma_nu`, each of length `length(y)`.
`mu_sigma` and `mu_nu` are exactly zero.

## See also

[`distrib_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
for the quantity this is the expectation of,
[`expected_hessian_exact.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.PseudoHuberDistrib.md)
for the predicate that reports this is not a closed form,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which reads that predicate, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)
eh <- distrib_expected_hessian(d, y, th)
vapply(eh, function(v) v[1], numeric(1))
#>        mu_mu  sigma_sigma        nu_nu     mu_sigma        mu_nu     sigma_nu 
#> -0.262082774 -0.917769722 -0.005412069  0.000000000  0.000000000 -0.066997583 

# The two entries odd in the residual are exactly zero by symmetry, so the
# location is orthogonal to the scale and the shape.
c(eh$mu_sigma[1], eh$mu_nu[1])
#> [1] 0 0

# Unlike the families that write their information out, this one reads
# `approx`: a Monte Carlo strategy gives a different, noisier answer. It is
# also the dear one, drawing from a generator that root-finds, so `nsim` is
# kept small here.
set.seed(1)
vapply(distrib_expected_hessian(d, y, th, approx = "mc", nsim = 200),
       function(v) v[1], numeric(1))
#>        mu_mu  sigma_sigma        nu_nu     mu_sigma        mu_nu     sigma_nu 
#> -0.283301188 -0.659597777 -0.001592754  0.000000000  0.000000000 -0.063362162 
```
