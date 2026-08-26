# Zero-Inflated Fourth Derivatives

The same two regimes as at third order, one order along. At \\y \> 0\\
the log-likelihood separates into \\\log(1-\zeta)\\ and the parent's, so
a component mixing \\\zeta\\ with a parent parameter is exactly zero; at
\\y = 0\\ the moment-to-cumulant sum over the partitions of a four-index
set gives \\d^I \log L_0\\, with every block naming \\\zeta\\ twice
contributing nothing because \\L_0\\ is affine in it.

The partitions of four indices number fifteen against five for three, so
this order is the more expensive of the two by that ratio and by nothing
else: the same identity and the same ratios \\d^B L_0 / L_0\\ are read.

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

  Logical of length 1. `TRUE` takes the expectation.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of fourth-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters and `zi`, each a numeric vector of length
`length(y)`. A one-parameter parent gives five.

## See also

[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
for the wrapper;
[`distrib_deriv3.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroInflatedDistrib.md)
for the order below;
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
the partition sum used here.

## Examples

``` r
zi <- zero_inflated(poisson_distrib())
th <- list(mu = 3, zi = 0.2)

names(distrib_deriv4(zi, 0, th))
#> [1] "mu_mu_mu_mu" "mu_mu_mu_zi" "mu_mu_zi_zi" "mu_zi_zi_zi" "zi_zi_zi_zi"
round(unlist(distrib_deriv4(zi, 0, th)), 4)
#> mu_mu_mu_mu mu_mu_mu_zi mu_mu_zi_zi mu_zi_zi_zi zi_zi_zi_zi 
#>      0.0234      0.1463      3.0823     81.5264  -1478.5073 

# Separation at a positive count holds at this order too.
distrib_deriv4(zi, 2, th)[["mu_mu_zi_zi"]]
#> [1] 0
```
