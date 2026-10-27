# The Score Identity Evaluated at Observations

Returns, for each component of the expected derivative of an order, the
quantity \\-\sum\_{\pi,\\\|\pi\| \ge 2} \prod\_{B \in \pi} \ell_B\\
evaluated at the observations `y` and the parameters `theta`, with no
expectation. It is the integrand of the score identity of
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md),
and it is what
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
reads and what
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
averages when `identity = "score"`.

## Usage

``` r
score_identity_terms(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 2 to 4.

## Value

A named list of component vectors, each of the length of `y`, named as
[`observed_deriv()`](https://statmodels7.github.io/distributions7/reference/observed_deriv.md)
names that order.

## See also

[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md),
[`expected_by_mc()`](https://statmodels7.github.io/distributions7/reference/expected_by_mc.md)
