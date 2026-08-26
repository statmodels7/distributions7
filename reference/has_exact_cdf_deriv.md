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
quadrature per component with two calls on the parent: measured on a
truncated gaussian's Hessian, 1.4 ms against 4.9 ms. Where the parent
has no closed form, though, that route differences its distribution
function and carries roughly \\10^{-8}\\ of relative error into the
Hessian, where the quadrature it replaced carried \\10^{-10}\\.

That is invisible in the Hessian itself, and visible downstream:
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)
differentiates the analytical Hessian, so a noisier Hessian degrades the
REFERENCE the fourth-order check compares against, and the check then
fails on code that is right.

## How a genuine method is recognized

`attr(m, "signature")[[1]]` is the class the method was registered on,
so an inherited fallback answers with a base class.
[`identical()`](https://rdrr.io/r/base/identical.html) on the method
object cannot be used for this, S7 wrapping it. The third- and
fourth-order defaults sit on `distrib` where the first two sit on
`continuous_distrib`, so both base classes are excluded or a stencil
would be read as a closed form.

The question is asked of the OWNING CLASS, so a family that registers a
method combining closed components with a stencil in one direction
answers `TRUE`. That is the intended reading: what the gate protects
against is the generic differencing the cdf itself.

## See also

[`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md),
its one caller, and
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the generics it asks about.

## Examples

``` r
# A discrete family answers TRUE at every order: its cdf derivatives are an
# exact sum.
vapply(1:4, function(k)
  distributions7:::has_exact_cdf_deriv(poisson_distrib(), k), logical(1))
#> [1] TRUE TRUE TRUE TRUE

# The gaussian writes all four out.
vapply(1:4, function(k)
  distributions7:::has_exact_cdf_deriv(gaussian1_distrib(), k), logical(1))
#> [1] TRUE TRUE TRUE TRUE

# The gamma does not, the derivative of an incomplete gamma in its shape
# having no elementary form, so truncation falls back to quadrature there.
vapply(1:4, function(k)
  distributions7:::has_exact_cdf_deriv(gamma2_distrib(), k), logical(1))
#> [1] FALSE FALSE FALSE FALSE
```
