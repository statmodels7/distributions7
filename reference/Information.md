# S7 Classes for the Information Matrix

`Information` is the abstract parent of the objects that select which
matrix is taken as the information of a likelihood: the observed
information, which is the negative Hessian of the log-likelihood at the
data, or the expected (Fisher) information. Its two subclasses are
`ObservedInfo` and `ExpectedInfo`. They are built by
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md),
which are the functions to call; the classes are exported so that
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html)
can test an object, and so that a package can define a subclass of
`Information` for a further choice of matrix.
[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
compares the two.

`ObservedInfo` has no properties. `ExpectedInfo` has one, `method`,
which is `NULL` or an object that inherits from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
and selects how the expected information is obtained. The raw
constructors do not validate their properties; the functions
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
do.

## Usage

``` r
Information()

ObservedInfo()

ExpectedInfo(method = NULL)
```

## Arguments

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
  for `ExpectedInfo`.

## Value

An S7 object of the class named, which inherits from `Information`.

## See also

[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
for the comparison of the two constructors;
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
for the constructors;
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)
for a function that takes the object.

## Examples

``` r
a <- expected_info(method = expect_mc(nsim = 500))
S7::S7_inherits(a, Information)
#> [1] TRUE
S7::S7_inherits(a, ExpectedInfo)
#> [1] TRUE
a@method
#> Expectation by Monte Carlo (nsim = 500)
```
