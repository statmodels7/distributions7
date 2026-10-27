# Expected Derivatives of a Multivariate Distribution, by Sampling

Computes the expected derivative of an order for a multivariate family
as an average over a sample drawn from the family. It serves the
fallback of
[`distrib_expected_hessian.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.multivariate_distrib.md)
and an explicit `method` given for a multivariate family that has a
closed form.

## Usage

``` r
mv_expected_sampled(distrib, y, theta, order, method)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- y:

  The observations; only the number of rows is used, because the
  expectation is over the law.

- theta:

  A named list of parameters, in which each component is a single
  number.

- order:

  The derivative order, 2 to 4.

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).

## Value

A named list of numeric vectors of length \\n\\, in which each vector is
constant.

## Details

A multivariate family has no quantiles, so
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
signals an error.
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
with `identity = "score"` average minus the products of the score
identity (at order 2, the negative outer product of the score) over
`method_nsim(method)` draws, and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
with `identity = "hessian"` averages the observed derivative of the
order. The averages agree in expectation and differ by sampling error at
any finite `nsim`; the sample size is 10000 for
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md).

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
