# Jacobian of the Inverse Link at the Estimate

Returns \\h_i'(\eta_i) = dg_i^{-1}/d\eta_i\\, one entry per parameter,
evaluated at the supplied linear predictors. Because each parameter
carries its own scalar link the Jacobian of \\\theta\\ in \\\eta\\ is
diagonal, and this is its diagonal.

It is what the delta method needs to carry a variance matrix from the
link scale, where the fit computes it, to the parameter scale, where it
is reported: \\\widehat{\mathrm{Var}}(\hat\theta) = J V J\\ with \\J =
\mathrm{diag}(h')\\. On a Gaussian at \\\eta = (1.5, \log 2.5)\\, with
the identity on \\\mu\\ and the logarithm on \\\sigma\\, it returns
`c(1, 2.5)`.

## Usage

``` r
fit_dtheta_deta(distrib, eta)
```

## Arguments

- distrib:

  An object inheriting from `distrib`, supplying `params` and
  `link_params`.

- eta:

  A numeric vector of linear predictors, one per parameter, in the order
  `distrib@params` gives. Only the first element of each link's answer
  is kept, so an `eta` longer than the parameter count is silently
  truncated by the same rule
  [`fit_eta_from_theta()`](https://statmodels7.github.io/distributions7/reference/fit_eta_from_theta.md)
  applies.

## Value

An unnamed numeric vector of length `length(distrib@params)`, in the
order `distrib@params` gives.

## See also

[`fit_theta_from_eta()`](https://statmodels7.github.io/distributions7/reference/fit_theta_from_eta.md),
the map this differentiates;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which uses it for the standard errors and nothing else;
[`linkfunctions7::dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.html).
