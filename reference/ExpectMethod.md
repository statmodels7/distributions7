# S7 Classes for the Method of an Expectation

`ExpectMethod` is the abstract parent of the objects that select how an
expected derivative is obtained. Its three subclasses are `ExpectOPG`,
`ExpectIntegrate` and `ExpectMC`. They are built by
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md),
which are the functions to call; the classes are exported so that
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html)
can test an object.
[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
compares the three.

`ExpectOPG` has no properties. `ExpectMC` has two, `nsim`, the number of
draws, and `identity`. `ExpectIntegrate` has seven, `identity` and the
tolerances and limits of the quadrature and of the series summation (see
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)).
The raw constructors `ExpectMC()` and `ExpectIntegrate()` do not
validate their properties, whereas
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
and
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
do.

## Usage

``` r
ExpectMethod()

ExpectOPG()

ExpectIntegrate(
  identity = "hessian",
  quad_atol = 1e-10,
  quad_rtol = 1e-08,
  max_panels = 4096L,
  series_atol = 1e-12,
  series_rtol = 1e-10,
  max_terms = 100000L
)

ExpectMC(nsim = integer(0), identity = "hessian")
```

## Arguments

- identity:

  One of `"hessian"` and `"score"`, the identity that the expectation is
  taken of, for `ExpectIntegrate` and `ExpectMC`. It is documented in
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- quad_atol, quad_rtol, max_panels, series_atol, series_rtol, max_terms:

  The properties of `ExpectIntegrate` that set the quadrature and the
  series summation, documented in
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md).

- nsim:

  The number of Monte Carlo draws, a single positive number, for
  `ExpectMC`.

## Value

An S7 object of the class named, which inherits from `ExpectMethod`.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
for the comparison of the three methods;
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md),
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
and
[`expect_mc()`](https://statmodels7.github.io/distributions7/reference/expect_mc.md)
for the constructors.

## Examples

``` r
a <- expect_mc(nsim = 500)
S7::S7_inherits(a, ExpectMethod)
#> [1] TRUE
S7::S7_inherits(a, ExpectMC)
#> [1] TRUE
a@nsim
#> [1] 500
a@identity
#> [1] "hessian"
```
