# The Link Scale Carried Back to Parameters

Applies each parameter's inverse link, \\\theta_i = g_i^{-1}(\eta_i)\\,
and returns the parameters as the named list every generic in the
package expects. It is the inverse of
[`fit_eta_from_theta()`](https://statmodels7.github.io/distributions7/reference/fit_eta_from_theta.md).

Every link maps onto the **open** interior of its parameter's domain,
and `linkfunctions7` clamps the result strictly inside when the
arithmetic saturates, so a \\\theta\\ obtained this way is admissible
whatever the optimizer proposed. On a Gaussian scale carried by the
logarithm, \\\eta = -800\\ gives \\\sigma = 1.9\times 10^{-77}\\ and
\\\eta = 800\\ gives \\1.8\times 10^{308}\\: both extreme, both still
inside \\(0, \infty)\\, so the density can be evaluated and the point
rejected on its likelihood. This is also why a confidence interval built
on the link scale and mapped back cannot run outside the domain.

## Usage

``` r
fit_theta_from_eta(distrib, eta)
```

## Arguments

- distrib:

  An object inheriting from `distrib`, supplying `params` and
  `link_params`.

- eta:

  A numeric vector of linear predictors, one per parameter, in the order
  `distrib@params` gives. A non-finite entry propagates to the parameter
  it names.

## Value

A named list of length `length(distrib@params)`, each component a number
on the parameter scale, named and ordered as `distrib@params`.

## See also

[`fit_eta_from_theta()`](https://statmodels7.github.io/distributions7/reference/fit_eta_from_theta.md),
the inverse;
[`fit_dtheta_deta()`](https://statmodels7.github.io/distributions7/reference/fit_dtheta_deta.md)
for the first derivative of the same map;
[`linkfunctions7::linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.html).
