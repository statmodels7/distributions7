# The Default Information of a Family

Returns the information object that a function takes when its
`information` argument is `NULL`:
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
when every parameter of the family has
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
of at least 2, that is, when the log-likelihood is twice differentiable
in every parameter, and
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
with `method = NULL` otherwise, for the whole matrix. A family with a
kink in a parameter (the Laplace, the elastic net) has a singular or
undefined observed information in that parameter, and the closed-form
expected information gives finite standard errors. A parameter whose
order is not established (`NA`) does not satisfy the condition. The rule
is the one of
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
and it is available to packages that take an `information` argument with
the same default.

## Usage

``` r
default_information(distrib)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

## Value

An
[Information](https://statmodels7.github.io/distributions7/reference/Information.md)
object: an `ObservedInfo` or an `ExpectedInfo` with `method = NULL`.

## See also

[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md),
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md),
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md),
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)

## Examples

``` r
default_information(gaussian1_distrib())
#> Observed information
default_information(laplace_distrib())
#> Expected information (closed form where the family has one, otherwise the outer product of the scores)
```
