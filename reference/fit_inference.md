# Standard Errors and Intervals From a Link-Scale Variance Matrix

Carries a variance matrix of the estimates on the link scale to the
quantities that a fit stores and reports: the variance matrix on the
parameter scale by the delta method, the standard errors on both scales
and the Wald intervals, which are symmetric on the link scale and mapped
through the inverse links.

## Usage

``` r
fit_inference(distrib, eta, V_eta, level)
```

## Arguments

- distrib:

  The fitted `distrib` object.

- eta:

  The estimates on the link scale.

- V_eta:

  The variance matrix on the link scale, named by the parameters;
  entries are `NA` where it could not be computed.

- level:

  The confidence level of the intervals.

## Value

A list with the components `vcov`, `vcov_eta`, `se`, `se_eta`, `ci` and
`ci_eta`, laid out as the properties of the same names in
[`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md).

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which stores the result for the observed information;
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
which recomputes it for another matrix.
