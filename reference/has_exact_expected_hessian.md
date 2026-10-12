# Does This Distribution Compute Its Expected Information Exactly?

`TRUE` when the distribution registers its own
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
method, rather than inheriting the one that approximates the
expectation.

## Usage

``` r
has_exact_expected_hessian(x)
```

## Arguments

- x:

  An object inheriting from class `"distrib"`.

## Value

A single logical.

## Details

The question decides whether the `approx` argument means anything. A
family with a closed form ignores it, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
rejects the argument in that case instead of accepting it and doing
something else: the Laplace is the example the package already
documents, where `approx` has no effect at all because
\\\mathcal{I}(\theta) = 1/b^{2}\\ is written out.

The comparison goes through
[`is_class()`](https://statmodels7.github.io/distributions7/reference/is_class.md)
rather than [`identical()`](https://rdrr.io/r/base/identical.html), for
the reason recorded there.

The argument is named `x` because the base class of the package is named
`distrib`, and an argument of that name would mask the class in the
comparison with the base class.

## See also

[`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md),
[`is_class()`](https://statmodels7.github.io/distributions7/reference/is_class.md)
