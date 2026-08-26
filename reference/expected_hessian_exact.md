# Is a Family's Expected Information Written Out?

A generic so that a family whose registered
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
method is not what its owning class suggests can say so;
[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md)
asks it and
[`expected_hessian_exact.distrib()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.distrib.md)
is the default.

The question
[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md)
asks, as a generic, so that a family whose registered method is not what
its owner suggests can say so.

## Usage

``` r
expected_hessian_exact(x, ...)
```

## Arguments

- x:

  An object inheriting from class `"distrib"`.

- ...:

  Passed to methods.

## Value

A single logical.

A single logical.

## Details

The default reads the class the
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
method is registered on: the base classes carry the approximating method
and every other registration is, by default, a family that wrote the
expectation out.

**Reading the owner is not sufficient, and two families prove it.** A
method registered on a family's own class may still be a CHAIN onto a
parent that approximates, and then the owner says "written out" about
arithmetic that is a quadrature. Measured at 100 observations, where
thirty-four families answer in a median of 0.183 ms:
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
costs 5220 ms, more than the
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
it chains onto, which costs 2230, and
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
costs 10980 ms. Both were reported as exact. The consequences were real
rather than cosmetic:
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
rejected a legitimate `fisher_scoring(approx = )` on those two with a
message stating that the family "computes its expected information in
closed form", which is untrue, and its standard-error branch entered a
multi-second quadrature believing it cheap.

A family that chains onto another therefore answers for its parent,
which is what the two methods registered here do.

## See also

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)

## Examples

``` r
distributions7:::expected_hessian_exact(gaussian1_distrib())
#> [1] TRUE
```
