# Register the Four CDF Derivative Orders of an Exponential Survival Family

Turns a function returning \\L = \log(1-F)\\ and its partial-derivative
evaluator into the four S7 methods, so that a family states its survival
function once instead of four times. Three families are registered
through it: the exponential, the Weibull and the generalized Pareto.

## Usage

``` r
register_surv_cdf(cls, pieces)
```

## Arguments

- cls:

  The S7 class to register on.

- pieces:

  A function of `(distrib, q, theta)` returning a list with `Lval` and
  `Lderiv`, and optionally `inside`, as
  [`surv_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/surv_cdf_deriv_k.md)
  documents.

## Value

Invisibly `NULL`. Called for the registration.

## Details

All four orders are registered,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
included, so these families take the closed route from the first order
up. Where the upper tail is asked for on the natural scale the
derivatives of \\S\\ are returned directly, \\\partial^I S = -\partial^I
F\\, which is the one case where the survival function is the quantity
the construction produces first.

`force(o)` inside the factory is what keeps the four registrations from
sharing one order.

## See also

[`surv_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/surv_cdf_deriv_k.md),
the body it registers;
[`distrib_grad_cdf.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.ExponentialDistrib.md),
[`distrib_grad_cdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Weibull1Distrib.md)
and
[`distrib_grad_cdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md),
the three families.
