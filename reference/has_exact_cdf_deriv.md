# Can the Parent Supply Exact CDF Derivatives?

`TRUE` when the parent has a genuine closed-form cdf derivative of the
given order, or is a discrete family whose cdf derivatives are an exact
sum.

## Usage

``` r
has_exact_cdf_deriv(parent, order)
```

## Arguments

- parent:

  The parent distribution.

- order:

  The cdf derivative order, 1 to 4.

## Value

A single logical.

## Details

Replacing the quadrature for \\d^B Z\\ with two calls on the parent's
cdf derivative is roughly an order of magnitude faster, and is taken
only where it is at least as accurate as what it replaces: when the
parent has no closed form that route differences its cdf, carrying more
relative error into the Hessian than the quadrature does.

That is invisible in the Hessian itself but not downstream:
[`numerical_deriv4`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)
differentiates the analytical Hessian, so a noisier Hessian degrades the
*reference* the fourth-order check compares against, and the check fails
on code that is right.

Detecting a genuine method uses the documented S7 trick:
`attr(m, "signature")[[1]]` is the class the method was registered on,
so an inherited fallback gives `continuous_distrib`.
[`identical()`](https://rdrr.io/r/base/identical.html) on the method
object does not answer the question, because S7 wraps it.

## See also

[`trunc_mass_derivs`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md)
