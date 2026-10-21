# Register the Third and Fourth Response Derivatives of a Family

Turns a function of `(distrib, y, theta, k)` into the two methods, so
that a family writes its rule once instead of twice.

## Usage

``` r
register_dy_k(cls, f)
```

## Arguments

- cls:

  The S7 class.

- f:

  The rule.

## Value

Invisibly `NULL`. It is called for the two registrations that it makes.

## Registered on

This body serves both third- and fourth-order response derivatives on
the families listed below, which write a rule in the response of their
own, so `?distrib_deriv3_y.Gamma1Distrib` and the other aliases open
this page: `Gamma1Distrib`, `Gamma2Distrib`, `ChisqDistrib`,
`ExponentialDistrib`, `Beta1Distrib`, `Beta2Distrib`, `Weibull1Distrib`,
`InvGauss1Distrib`, `InvGauss2Distrib`, `Lognormal1Distrib`,
`VonMises1Distrib`, `VonMises2Distrib`.

## See also

[`loc_deriv_y_k()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md),
the companion for the location families;
[`distrib_deriv3_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.continuous_distrib.md)
for the stencil to which a family outside both lists falls back;
[`dy_log()`](https://statmodels7.github.io/distributions7/reference/dy_log.md)
for the elementary terms from which each rule is built.
