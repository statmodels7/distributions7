# The Kink a Family Carries

The one non-smooth piece of a log-density, as the composition
\\c(\theta)\\\phi(v(y,\theta))\\ of
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md),
or `NULL` for a family that is smooth in every parameter.
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
reads the order of differentiability off the answer.

## Usage

``` r
kink_decomposition(distrib)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

## Value

A
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md),
or `NULL` when the family is smooth.

## Details

Three families declare one, all for the same reason – an absolute value
of the residual in the log-density – and every other family inherits the
base method, which returns `NULL`.

A wrapper does not propagate the declaration today, so a wrapper of a
kinked family returns `NULL` here while `params_smooth` still records
the kink;
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
reports `NA` there rather than `Inf`, which says the order has not been
established rather than that the parameter is smooth.

## See also

[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md)
for the object,
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
for the order deduced from it,
[`check_kink()`](https://statmodels7.github.io/distributions7/reference/check_kink.md)
for the validator, and
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md)
for the logical this generalises.

## Examples

``` r
kink_decomposition(gaussian1_distrib())
#> NULL
kink_decomposition(laplace_distrib())
#> <kink_spec> c(theta) * phi(v),  phi = abs,  order 0
```
