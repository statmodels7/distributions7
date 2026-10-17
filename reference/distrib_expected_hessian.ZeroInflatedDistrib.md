# Zero-Inflated Expected Information

Computes the expectation of the observed Hessian in closed form, by
splitting the expectation over the two events \\y = 0\\ and \\y \> 0\\.
The zero contributes with probability \\L_0\\ and carries the mixture's
corrections; the positive part contributes with probability \\1 - L_0\\
and carries the parent's own expected Hessian, less what the zero would
have contributed to it. The wrapper's own blocks do not depend on the
data, so the returned vectors are constant when the parent's expected
Hessian is, and `y` is then used only for its length.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- y:

  A numeric vector of observations. Only its length is used when the
  parent's own expected information is closed form; where it is not, the
  parent block reads `y` itself, through
  `distrib_expected_hessian(parent, y, ...)`.

- theta:

  A named list with the parent's parameters followed by `zi`.

- scale:

  One of `"parameter"` (the default) or `"link"`. The generic applies
  the link scale to the value that the method returns.

- approx:

  Forwarded to the parent's
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
  and read only where the parent has no closed form for it; every other
  block is an exact combination of the parent's density, score and
  observed Hessian at zero.
  [`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
  returns, for this class, the value that it returns for the parent.

- nsim:

  Forwarded for the same reason and under the same condition.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length `length(y)`, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
Constant across observations when the parent's own expected information
is; a parent evaluated by an approximation can vary by observation.

## Details

The blocks of the wrapper are exact combinations of the parent's mass,
score and observed Hessian at zero and of the parent's expected Hessian,
an expectation over the whole support from which the single point at
zero is subtracted in closed form. `approx` and `nsim` are forwarded to
the parent's
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
and take effect only where the parent has no closed form, so the result
is exact exactly when the parent's is, as
[`expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
reports.

## Notation

\\f\\ is the parent's mass function, \\\zeta\\ the probability of a
structural zero, \\L_0 = \zeta + (1-\zeta)f(0)\\ the inflated mass at
zero, \\w = (1-\zeta)f(0)/L_0\\ the posterior probability that an
observed zero came from the parent, \\s\\ the parent's score, \\H\\ its
observed Hessian and \\\ell\\ the log-mass of one observation.

## See also

[`distrib_hessian.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroInflatedDistrib.md)
for the observed matrix,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
whose Fisher scoring inverts this, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)
set.seed(2)
y <- distrib_rng(d, 200, theta)

EH <- distrib_expected_hessian(d, y, theta)
vapply(EH, function(z) z[1], numeric(1))
#>      mu_mu      zi_zi      mu_zi 
#> -0.2175121 -4.4092338  0.1732687 

# Closed form, so two calls agree to the bit and the method draws no samples.
identical(EH, distrib_expected_hessian(d, y, theta))
#> [1] TRUE

# The closed form equals the observed Hessian summed against the mass
# function, over a support that extends far enough out.
sup <- 0:400
m <- distrib_pdf(d, sup, theta)
Hs <- distrib_hessian(d, sup, theta)
rbind(summed = vapply(Hs, function(z) sum(z * m), numeric(1)),
      closed = vapply(EH, function(z) z[1], numeric(1)))
#>             mu_mu     zi_zi     mu_zi
#> summed -0.2175121 -4.409234 0.1732687
#> closed -0.2175121 -4.409234 0.1732687
```
