# Is a Family's Expected Information Written Out?

`TRUE` when the family computes its expected information in closed form,
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

Six of the shipped univariate families answer `FALSE`:
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md),
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md),
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
and
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

Asking the OWNING CLASS of the registered method is not enough on its
own, which is why the generic exists: the pseudo-Huber registers a
method of its own that calls the fallback and then patches the two
components vanishing by symmetry, and the centered skew normal chains
onto a parent whose expected information is itself a quadrature. Both
would read as exact from the method's owner alone, and both cost seconds
where a closed form costs milliseconds.

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
#> [1] FALSE
```
