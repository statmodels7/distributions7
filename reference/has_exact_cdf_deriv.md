# Can the Parent Supply Exact CDF Derivatives?

Reports whether the parent has a genuine closed-form derivative of its
distribution function at the given order, or is a discrete family, whose
cdf derivatives are an exact finite sum.
[`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md)
takes the cheap route only where this is `TRUE`, and falls back to
quadrature otherwise.

## Usage

``` r
has_exact_cdf_deriv(parent, order)
```

## Arguments

- parent:

  The parent distribution, before wrapping.

- order:

  The derivative order, an integer from 1 to 4.

## Value

A single logical.

## Why the cheap route is gated

Reading \\d^B Z\\ off the parent's cdf derivatives replaces one
quadrature per component with two calls on the parent. Where the parent
has no closed-form cdf derivative, that route differences its
distribution function and has a larger relative error than the
quadrature that it would replace. The error passes into
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md),
which differentiates the analytical Hessian, and enlarges the error of
the reference used in the fourth-order check.

## How a genuine method is recognized

`attr(m, "signature")[[1]]` is the class on which the method was
registered, so for an inherited fallback it is a base class.
[`identical()`](https://rdrr.io/r/base/identical.html) on the method
object cannot be used for this, because S7 wraps the method. The third-
and fourth-order defaults sit on `distrib` where the first two sit on
`continuous_distrib`, so both base classes are excluded or a stencil
would be read as a closed form.

The test uses the class that owns the method, so for a family that
registers a method combining closed components with a stencil in one
direction the result is `TRUE`. This behavior is intended, because the
gate protects against the generic method differencing the cdf itself.

## See also

[`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md),
its caller in this file, and
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the generics that it tests.

## Examples

``` r
# For a discrete family the result is TRUE at every order, because its cdf
# derivatives are an exact sum.
vapply(1:4, function(k)
  distributions7:::has_exact_cdf_deriv(poisson_distrib(), k), logical(1))
#> [1] TRUE TRUE TRUE TRUE

# The gaussian has closed forms for all four orders.
vapply(1:4, function(k)
  distributions7:::has_exact_cdf_deriv(gaussian1_distrib(), k), logical(1))
#> [1] TRUE TRUE TRUE TRUE

# The gamma's cdf derivatives are compiled, so the result is TRUE as well.
vapply(1:4, function(k)
  distributions7:::has_exact_cdf_deriv(gamma2_distrib(), k), logical(1))
#> [1] TRUE TRUE TRUE TRUE
```
