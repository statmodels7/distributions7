# Zero-Adjusted Third Derivatives, Continuous Parent

A continuous parent puts **no** mass at zero, so there is nothing to
truncate away and no normalizing constant to differentiate. Away from
the atom the log-likelihood is \\\log(1 - \pi) + \ell(y;\theta)\\: the
\\\theta\\ derivatives are the parent's unchanged, and every component
mixing `za` with a parent parameter is exactly zero.

The object is a **mixed** distribution, a density on the positive line
plus an atom at zero, and it declares that atom through
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md),
so
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
tests it correctly: the density integrates to \\1 - \pi\\, not to 1.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations, zero included.

- theta:

  A named list with the parent's parameters followed by `za`, the
  probability of the atom, in \\(0, 1)\\.

- expected:

  Logical of length 1. `TRUE` takes the expectation over the mixed law,
  atom included.

- ...:

  Passed on, `approx` and `nsim` among them.

## Value

A named list of third-derivative components, keyed lexicographically by
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
on the parent's parameters and `za`, each a numeric vector of length
`length(y)`. A two-parameter parent gives ten.

## See also

[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the wrapper;
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md),
which declares the atom;
[`distrib_deriv4.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ZeroAdjustedContinuousDistrib.md)
for the order above;
[`distrib_deriv3.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.ZeroAdjustedDiscreteDistrib.md),
where a truncation constant does appear.

## Examples

``` r
za <- zero_adjusted(gamma2_distrib())
th <- list(mu = 2, sigma2 = 0.7, za = 0.25)

names(distrib_deriv3(za, c(0, 1.5, 3), th))
#>  [1] "mu_mu_mu"             "mu_mu_sigma2"         "mu_mu_za"            
#>  [4] "mu_sigma2_sigma2"     "mu_sigma2_za"         "mu_za_za"            
#>  [7] "sigma2_sigma2_sigma2" "sigma2_sigma2_za"     "sigma2_za_za"        
#> [10] "za_za_za"            

# Away from the atom the theta derivatives are the parent's, exactly.
all.equal(distrib_deriv3(za, c(1.5, 3), th)[["mu_mu_mu"]],
          distrib_deriv3(gamma2_distrib(), c(1.5, 3),
                         list(mu = 2, sigma2 = 0.7))[["mu_mu_mu"]])
#> [1] TRUE

# And the mixed components vanish.
distrib_deriv3(za, c(0, 1.5, 3), th)[["mu_mu_za"]]
#> [1] 0 0 0
```
