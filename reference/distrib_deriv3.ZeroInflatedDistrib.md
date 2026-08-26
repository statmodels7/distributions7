# Zero-Inflated Third Derivatives

Two regimes, and the observation decides which. At \\y \> 0\\ the
log-likelihood is \\\log(1-\zeta) + \ell(y;\theta)\\, which
**separates**: every mixed component in \\\zeta\\ and a parent parameter
is exactly zero. At \\y = 0\\ it is \\\log L_0\\ with \\L_0 = \zeta +
(1-\zeta) f(0;\theta)\\, and the derivatives follow from the
moment-to-cumulant expansion over set partitions, \$\$d^I \log L_0 =
\sum\_{\pi} (-1)^{\|\pi\|-1}(\|\pi\|-1)! \prod\_{B \in \pi} \frac{d^B
L_0}{L_0}.\$\$

\\L_0\\ is **affine** in \\\zeta\\, so any block of a partition naming
\\\zeta\\ twice contributes nothing and the sum is shorter than it
looks. A whole component may still be non-zero: its partition into
singletons survives.

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- y:

  A numeric vector of counts.

- theta:

  A named list with the parent's parameters followed by `zi`, the
  inflation probability in \\(0, 1)\\.

- expected:

  Logical of length 1. `TRUE` takes the expectation, which for a lattice
  parent is an exact sum over the support.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of third-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters and `zi`, each a numeric vector of length
`length(y)`.

## See also

[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
for the wrapper and for why it cannot be stacked;
[`distrib_deriv4.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ZeroInflatedDistrib.md)
for the order above;
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
the partition sum used here.

## Examples

``` r
zi <- zero_inflated(poisson_distrib())
th <- list(mu = 3, zi = 0.2)

# At a positive count the likelihood separates, so a mixed component is 0.
distrib_deriv3(zi, 2, th)[["mu_zi_zi"]]
#> [1] 0

# At zero it is not: log L0 mixes the two.
round(unlist(distrib_deriv3(zi, 0, th)), 4)
#> mu_mu_mu mu_mu_zi mu_zi_zi zi_zi_zi 
#>  -0.0925  -0.5781  -6.8590 124.3896 
```
