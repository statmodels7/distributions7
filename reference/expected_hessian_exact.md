# Is a Family's Expected Information Written Out?

`TRUE` when the family computes its expected information exactly, in
closed form or by a quadrature of its own to working precision, and
`FALSE` when a call to
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
reaches a fallback and the answer is therefore an approximation. It is a
generic so that a family whose registered method is not what its owning
class suggests can say so;
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

The predicate is what a consumer reads to decide which information to
report. Where it answers `TRUE` the expected information costs one
evaluation and has the smaller variance, so it is the better matrix to
invert; where it answers `FALSE` the choice is between an approximation
and the observed Hessian, which every family has and which is exact.
That is the rule
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
follows for its standard errors and the one statmodels7 follows in
[`vcov()`](https://rdrr.io/r/stats/vcov.html).

Two of the shipped univariate families answer `FALSE`:
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
and
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).
The skew normals, the skew t and the pseudo-Huber answer `TRUE`, their
expected information being one quadrature over the standardized response
per distinct shape; see
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md).

Asking the OWNING CLASS of the registered method is not enough on its
own, which is why the generic exists: a method registered on a family's
own class may chain onto a parent that approximates, and a family that
does so answers for its parent.

The default reads the class the
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
method is registered on: the base classes carry the approximating method
and every other registration is, by default, a family that wrote the
expectation out.

**Reading the owner is not sufficient in general.** A method registered
on a family's own class may be a CHAIN onto a parent that approximates,
and then the owner says "written out" about arithmetic that is not. A
family that chains onto another therefore answers for its parent, which
is what the wrappers and
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
do, the last answering for
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

## See also

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
the wrapper that asks this;
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the quantity itself, and
[`expected_by_opg()`](https://statmodels7.github.io/distributions7/reference/expected_by_opg.md)
for the fallback a `FALSE` answer reaches.

[`has_exact_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/has_exact_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)

## Examples

``` r
# A Gaussian writes its information out; a Poisson-inverse gaussian does not.
expected_hessian_exact(gaussian1_distrib())
#> [1] TRUE
expected_hessian_exact(pig1_distrib())
#> [1] TRUE
```
