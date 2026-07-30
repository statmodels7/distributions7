# Map the Link Scale Back to Parameters

Applies each parameter's inverse link, \\\theta_i = g_i^{-1}(\eta_i)\\.

## Usage

``` r
fit_theta_from_eta(distrib, eta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- eta:

  A numeric vector of linear predictors, one per parameter.

## Value

A named list of parameters on the natural scale.

## Details

The inverse of
[`fit_eta_from_theta`](https://statmodels7.github.io/distributions7/reference/fit_eta_from_theta.md).
Because every link maps onto its parameter's domain by construction, a
\\\theta\\ obtained this way is admissible whatever the optimizer
proposed – which is the reason for working on the link scale at all, and
the reason a confidence interval built there and mapped back cannot run
outside the domain.

## See also

[`fit_eta_from_theta`](https://statmodels7.github.io/distributions7/reference/fit_eta_from_theta.md)
