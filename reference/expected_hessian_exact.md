# Whether a Family's Expected Information Is Written Out

Returns `TRUE` when the family computes its expected information
exactly, in closed form or by a quadrature of its own to working
precision, and `FALSE` when a call to
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
reaches a fallback and the answer is therefore an approximation. It is a
generic so that a family whose registered method is not what its owning
class suggests can report this;
[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md)
calls it and
[`expected_hessian_exact.distrib()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.distrib.md)
is the default.

Implements the test that
[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md)
applies, as a generic so that a family whose registered method differs
from what its owning class indicates can supply its own answer.

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

A function that reports the information of a family uses the predicate
to choose which information to report. Where it returns `TRUE` the
expected information costs one evaluation and has the smaller variance,
so it is the better matrix to invert; where it returns `FALSE` the
choice is between an approximation and the observed Hessian, which every
family has and which is exact.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
follows this rule for its standard errors, and statmodels7 follows it in
[`vcov()`](https://rdrr.io/r/stats/vcov.html).

Every shipped univariate family returns `TRUE`. The Poisson-inverse
Gaussian families sum their expected information exactly over the
support, and the skew normals, the skew t and the pseudo-Huber evaluate
one quadrature over the standardized response per distinct shape; see
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md).
Among the wrappers,
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
returns `FALSE`.

Reading the class that owns the registered method is not sufficient on
its own, and the generic exists for this reason: a method registered on
the class of a family may chain onto a parent that approximates, and a
family that does so returns the answer of its parent.

The default reads the class on which the
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
method is registered: the base classes carry the approximating method
and every other registration is, by default, a family that wrote the
expectation out.

Reading the owner is not sufficient in general. A method registered on
the class of a family may be a chain onto a parent that approximates,
and then the owner reads as "written out" although the arithmetic is
not. A family that chains onto another therefore returns the answer of
its parent, as the wrappers do;
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
returns `FALSE` unconditionally.

## See also

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
the wrapper that calls it;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the quantity itself, and
[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md)
for the fallback that a `FALSE` answer reaches.

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)

## Examples

``` r
# A Gaussian writes its information out; a truncated Gaussian does not.
expected_hessian_exact(gaussian1_distrib())
#> [1] TRUE
expected_hessian_exact(truncated(gaussian1_distrib(), lower = 0))
#> [1] FALSE
```
