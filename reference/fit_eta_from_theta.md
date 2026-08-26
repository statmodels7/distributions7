# Parameters Carried to the Link Scale

Applies each parameter's own link, \\\eta_i = g_i(\theta_i)\\, and
returns the resulting vector of linear predictors. This is how a
starting value expressed in natural parameters enters the optimizer,
which works on \\\eta\\ because it is unconstrained.

Only the first element of each component of `theta` is read. A fit
estimates one \\\theta\\ for the whole sample, so a component of length
\\n\\ is a vector the caller built for a density evaluation, not \\n\\
separate parameters: `list(mu = c(1, 99, 99), sigma = 2)` gives the same
\\\eta\\ as `list(mu = 1, sigma = 2)`.

## Usage

``` r
fit_eta_from_theta(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from `distrib`, supplying `params` and
  `link_params`.

- theta:

  A named list of parameters on the parameter scale, in any order; the
  components are read by position in `distrib@params`, so it must
  already have been through
  [`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md).
  A value outside the parameter's domain reaches the link, which returns
  `NaN` or a non-finite value there.

## Value

An unnamed numeric vector of length `length(distrib@params)`, one linear
predictor per parameter, in the order `distrib@params` gives.

## See also

[`fit_theta_from_eta()`](https://statmodels7.github.io/distributions7/reference/fit_theta_from_eta.md),
the inverse;
[`fit_dtheta_deta()`](https://statmodels7.github.io/distributions7/reference/fit_dtheta_deta.md)
for the Jacobian between the two scales;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
the caller.
